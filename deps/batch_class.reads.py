#!/usr/bin/env python3
"""
宏基因组/病毒组终极统一分类框架 (Metagenome Batch Classifier) - 工业级美学版
特色: 
  - [参数支持] 支持 --logs-dir 分离日志，支持手动指派每个工具的数据库
  - [输出架构] 样本级绝对压平化 (Flatten) 输出，告别多层子目录
  - [工业 UI ] 屏幕输出极度优美，带任务编号追踪与底部进度条锁定
  - [智能审计] 自动捕获时间、峰值内存，断点续传时无损继承历史资源记录
"""

import os
import sys
import argparse
import logging
import subprocess
import time
import re
import shutil
import csv
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# ==========================================
# 1. 工业级 UI 与日志系统
# ==========================================
class TqdmLoggingHandler(logging.Handler):
    """确保日志打印不会冲断底部的 tqdm 进度条"""
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

def setup_logging(logs_dir):
    logs_path = Path(logs_dir)
    pipeline_log_dir = logs_path / "pipeline_logs"
    pipeline_log_dir.mkdir(parents=True, exist_ok=True)
    log_file = pipeline_log_dir / "master_pipeline.log"
    
    logger = logging.getLogger("Master")
    logger.setLevel(logging.INFO)
    logger.handlers = [] # 清除默认
    
    # 后台文件日志
    fh = logging.FileHandler(log_file, mode='a')
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    
    # 屏幕高亮日志
    ch = TqdmLoggingHandler(level=logging.INFO)
    ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger, log_file

def find_samples(input_dir):
    input_path = Path(input_dir)
    extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
    all_files = []
    for ext in extensions:
        all_files.extend(list(input_path.glob(f'*{ext}')))
        
    samples = {}
    processed = set()
    paired_pattern = re.compile(r'(.*)([_.]R?)(1)([_.]?(?:001)?(?:_unmapped)?\.(?:fastq|fq)(?:\.gz)?)$', re.IGNORECASE)
    
    for file_path in all_files:
        if file_path in processed: continue
        match = paired_pattern.match(file_path.name)
        if match:
            base_name, sep, _, suffix = match.groups()
            r2_name = f"{base_name}{sep}2{suffix}"
            r2_path = input_path / r2_name
            if r2_path in all_files and r2_path not in processed:
                samples[base_name] = {'type': 'paired', 'r1': str(file_path), 'r2': str(r2_path)}
                processed.update([file_path, r2_path])
                continue
                
    for file_path in all_files:
        if file_path not in processed:
            sname = file_path.name
            for suf in extensions + ['.unmapped']: sname = sname.replace(suf, '')
            sname = re.sub(r'[_.]R?1$', '', sname)
            samples[sname] = {'type': 'single', 'r1': str(file_path), 'r2': None}
            processed.add(file_path)
            
    return samples

# ==========================================
# 2. 数据库智能极速发现 (支持 .cfr 与手动指派)
# ==========================================
def resolve_databases(args, req_tools, logger):
    dbs = {}
    base_dir = Path(args.db_dir) if args.db_dir else None
    
    manual_args = {
        'kraken2': args.db_kraken2, 'kraken2x': args.db_kraken2x, 'krakenuniq': args.db_krakenuniq,
        'centrifuger': args.db_centrifuger, 'ganon': args.db_ganon, 'kaiju': args.db_kaiju,
        'kunpeng': args.db_kunpeng, 'metabuli': args.db_metabuli, 'sylph': args.db_sylph
    }

    def extract_prefix_or_file(tool, path_obj):
        if not path_obj.is_dir(): return str(path_obj)
        if tool == 'centrifuger':
            cfs = list(path_obj.glob("*.1.cfr")) # 严格匹配 .cfr
            return str(cfs[0]).replace(".1.cfr", "") if cfs else str(path_obj)
        elif tool == 'ganon':
            ibfs = list(path_obj.glob("*.ibf"))
            return str(ibfs[0]).replace(".ibf", "") if ibfs else str(path_obj)
        elif tool == 'kaiju':
            fmis = list(path_obj.glob("*.fmi"))
            return str(fmis[0]) if fmis else str(path_obj)
        elif tool == 'sylph':
            syldbs = list(path_obj.glob("*.syldb"))
            # 选 c 值最大的数据库（c越大越兼容高深度reads）
            if syldbs:
                import re
                def _c_val(p):
                    m = re.search(r'\.c(\d+)\.', p.name)
                    return int(m.group(1)) if m else 0
                syldbs.sort(key=_c_val, reverse=True)
                return str(syldbs[0])
            return str(path_obj)
        return str(path_obj)

    for tool in req_tools:
        # 1. 优先使用手动指定的路径
        manual_path = manual_args.get(tool)
        if manual_path and Path(manual_path).exists():
            dbs[tool] = extract_prefix_or_file(tool, Path(manual_path))
            continue
            
        # 2. 浅层极速搜索 --db-dir
        if base_dir and base_dir.exists():
            for item in base_dir.iterdir():
                if not item.is_dir(): continue
                item_name = item.name.lower()
                
                if tool == 'kraken2' and ('kraken2' in item_name and 'x' not in item_name):
                    dbs[tool] = extract_prefix_or_file(tool, item)
                    break
                elif tool == 'kraken2x' and ('kraken2-x' in item_name or 'kraken2x' in item_name):
                    dbs[tool] = extract_prefix_or_file(tool, item)
                    break
                elif tool != 'kraken2' and tool != 'kraken2x':
                    search_key = tool.replace('_', '')
                    item_compact = item_name.replace('_', '')
                    # kaiju 特殊处理：移除 db 后缀干扰
                    item_compact = item_compact.replace('kaijudb', 'kaiju')
                    if search_key in item_compact:
                        dbs[tool] = extract_prefix_or_file(tool, item)
                        break

    if 'kaiju' in dbs:
        if args.kaiju_nodes and Path(args.kaiju_nodes).exists():
            dbs['kaiju_nodes'] = args.kaiju_nodes
            dbs['kaiju_names'] = args.kaiju_names
        else:
            k_parent = Path(dbs['kaiju']).parent
            search_dirs = [k_parent, k_parent.parent]
            if base_dir: search_dirs.append(base_dir / "taxonomy")
            for d in search_dirs:
                if (d / "nodes.dmp").exists() and (d / "names.dmp").exists():
                    dbs['kaiju_nodes'] = str(d / "nodes.dmp")
                    dbs['kaiju_names'] = str(d / "names.dmp")
                    break

    return dbs

# ==========================================
# 3. 架构基石：输出压平与隔离日志
# ==========================================
class BaseTool:
    def __init__(self, name, outdir_base, logs_dir_base, threads, db_path):
        self.name = name
        self.outdir_base = Path(outdir_base)
        self.logs_dir_base = Path(logs_dir_base)
        self.threads = threads
        self.db_path = db_path

    def _get_sample_dirs(self, sname):
        s_dir = self.outdir_base / sname
        log_dir = self.logs_dir_base / sname
        s_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        return s_dir, log_dir

    def run_cmd(self, cmd, log_file, redirect_stdout=None):
        start_time = time.time()
        time_bin = shutil.which("time")
        is_linux = sys.platform.startswith('linux')
        
        if time_bin and is_linux: 
            cmd = [time_bin, "-v"] + cmd
            
        success = False
        try:
            with open(log_file, 'w') as f_log:
                f_log.write(f"=== [{self.name.upper()}] COMMAND ===\n{' '.join(cmd)}\n{'='*40}\n\n")
                f_log.flush()
                
                if redirect_stdout:
                    with open(redirect_stdout, 'w') as f_out:
                        process = subprocess.run(cmd, stdout=f_out, stderr=f_log, text=True)
                else:
                    process = subprocess.run(cmd, stdout=f_log, stderr=subprocess.STDOUT, text=True)
                    
                success = (process.returncode == 0)
                if not success:
                    f_log.write(f"\n\n[ERROR] Command failed with exit code: {process.returncode}\n")
        except Exception as e:
            with open(log_file, 'a') as f_log: 
                f_log.write(f"\n\n[PYTHON EXCEPTION]\n{traceback.format_exc()}\n")
            success = False
            
        elapsed_time = time.time() - start_time
        max_rss_mb = 0.0
        cpu_pct = 0.0
        
        if os.path.exists(log_file) and is_linux:
            try:
                with open(log_file, 'r') as f:
                    txt = f.read()
                    matches = re.findall(r'Maximum resident set size \(kbytes\):\s+(\d+)', txt)
                    if matches: max_rss_mb = max([float(m) for m in matches]) / 1024.0
                    m_cpu = re.search(r'Percent of CPU this job got:\s+(\d+)', txt)
                    if m_cpu: cpu_pct = float(m_cpu.group(1))
            except: pass
                
        return success, elapsed_time, max_rss_mb, cpu_pct

    def safe_process(self, sname, info):
        try:
            return self.process_sample(sname, info)
        except Exception as e:
            _, log_dir = self._get_sample_dirs(sname)
            with open(log_dir / f"{self.name}.CRASH.log", 'w') as f:
                f.write(traceback.format_exc())
            return sname, False, 0.0, 0.0, 0.0, "Crash!"

    def process_sample(self, sname, info):
        raise NotImplementedError

# --- 工具子类断点续传与压平化逻辑 ---
class Kraken2Tool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        out_res = s_dir / f"{sname}.kraken2.out"
        out_rep = s_dir / f"{sname}.kraken2.report"
        if out_rep.exists() and out_rep.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "kraken2 ..."
        
        cmd = ["kraken2", "--db", self.db_path, "--report", str(out_rep), "--output", str(out_res), 
               "--threads", str(self.threads), "--report-minimizer-data"]
        if info['type'] == 'paired': cmd.extend(["--paired", info['r1'], info['r2']])
        else: cmd.append(info['r1'])
        
        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log")
        if not succ and out_res.exists(): out_res.unlink()
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

class Kraken2XTool(BaseTool):
    def _extract_report(self, min_rep, final_rep):
        with open(min_rep, 'r') as fin, open(final_rep, 'w') as fout:
            for line in fin:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 8: fout.write(f"{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t{c[6]}\t{c[7]}\n")

    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        k2x_out = s_dir / f"{sname}.kraken2-x.out"
        min_rep = s_dir / f"{sname}.kraken2-x.minimizer.report"
        final_rep = s_dir / f"{sname}.kraken2-x.report"
        if final_rep.exists() and final_rep.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "kraken2 (x) ..."
        
        cmd = ["kraken2", "--db", self.db_path, "--threads", str(self.threads),
               "--report", str(min_rep), "--output", str(k2x_out), "--report-minimizer-data"]
        if info['type'] == 'paired': cmd.extend(["--paired", info['r1'], info['r2']])
        else: cmd.append(info['r1'])
        
        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log")
        if succ and min_rep.exists():
            start_ext = time.time()
            self._extract_report(min_rep, final_rep)
            t += (time.time() - start_ext)
        else:
            for f in [k2x_out, min_rep, final_rep]:
                if f.exists(): f.unlink()
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

class KrakenUniqTool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        out_res = s_dir / f"{sname}.krakenuniq.out"
        out_rep = s_dir / f"{sname}.krakenuniq.report.txt"
        if out_rep.exists() and out_rep.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "krakenuniq ..."
        
        cmd = ["krakenuniq", "--db", self.db_path, "--threads", str(self.threads),
               "--output", str(out_res), "--report-file", str(out_rep), "--preload"]
        if info['type'] == 'paired': cmd.extend(["--paired", info['r1'], info['r2']])
        else: cmd.append(info['r1'])
        
        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log")
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

class CentrifugerTool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        out_class = s_dir / f"{sname}.centrifuger.out"
        out_quant = s_dir / f"{sname}.centrifuger.quant.tsv"
        log_f = log_dir / f"{self.name}.log"
        if out_quant.exists() and out_quant.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "centrifuger ..."
        
        cmd1 = ["centrifuger", "-x", self.db_path, "-t", str(self.threads)]
        if info['type'] == 'paired': cmd1.extend(["-1", info['r1'], "-2", info['r2']])
        else: cmd1.extend(["-u", info['r1']])
        
        s1, t1, m1, _c1 = self.run_cmd(cmd1, log_f, redirect_stdout=out_class)
        if not s1: 
            if out_class.exists(): out_class.unlink()
            return sname, False, t1, m1, 0.0, ' '.join(cmd1)
        
        cmd2 = ["centrifuger-quant", "-x", self.db_path, "-c", str(out_class)]
        s2, t2, m2, _c2 = self.run_cmd(cmd2, str(log_f)+".quant", redirect_stdout=out_quant)
        
        with open(log_f, 'a') as f1, open(str(log_f)+".quant", 'r') as f2:
            f1.write("\n\n" + f2.read())
        Path(str(log_f)+".quant").unlink()

        if not s2 and out_quant.exists(): out_quant.unlink()
        return sname, s2, (t1+t2), max(m1, m2), 0.0, f"{' '.join(cmd1)} && centrifuger-quant..."

class GanonTool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        prefix = s_dir / f"{sname}.ganon"
        rep = Path(f"{prefix}.rep")
        if rep.exists() and rep.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "ganon ..."
        
        cmd = ["ganon", "classify", "--db-prefix", self.db_path, "--output-prefix", str(prefix), "--threads", str(self.threads)]
        if info['type'] == 'paired': cmd.extend(["--paired-reads", info['r1'], info['r2']])
        else: cmd.extend(["--single-reads", info['r1']])
        
        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log")
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

class KunpengTool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        succ_mark = s_dir / f"{sname}.kunpeng.success"
        if succ_mark.exists(): return sname, True, 0.0, 0.0, 0.0, "kunpeng ..."
        
        # 使用 tmp 文件夹防止污染，成功后压平移动
        kp_tmp = s_dir / f"{sname}_kunpeng_tmp"
        if kp_tmp.exists(): shutil.rmtree(kp_tmp)
        kp_tmp.mkdir(parents=True, exist_ok=True)
        
        cmd = ["kun_peng", "direct", "--db", self.db_path, "--num-threads", str(self.threads), "--output-dir", str(kp_tmp)]
        if info['type'] == 'paired': cmd.extend(["--paired-end-processing", info['r1'], info['r2']])
        else: cmd.append(info['r1'])
        
        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log")
        
        if succ:
            for f in kp_tmp.iterdir():
                shutil.move(str(f), str(s_dir / f"{sname}.kunpeng.{f.name}"))
            shutil.rmtree(kp_tmp)
            succ_mark.touch()
        else: 
            shutil.rmtree(kp_tmp)
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

class KaijuTool(BaseTool):
    def __init__(self, outdir, logs_dir, threads, dbs):
        super().__init__('kaiju', outdir, logs_dir, threads, dbs.get('kaiju'))
        self.nodes = dbs.get('kaiju_nodes')
        self.names = dbs.get('kaiju_names')
        
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        out_raw = s_dir / f"{sname}.kaiju.out"
        out_names = s_dir / f"{sname}.kaiju.names.out"
        out_sum = s_dir / f"{sname}.kaiju.summary.tsv"
        log_f = log_dir / f"{self.name}.log"
        
        if out_sum.exists() and out_sum.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "kaiju ..."
        if not self.nodes or not self.names: 
            with open(log_f, 'w') as f: f.write("Error: Missing nodes.dmp or names.dmp\n")
            return sname, False, 0.0, 0.0, 0.0, "kaiju (missing db)"
        
        cmd1 = ["kaiju", "-z", str(self.threads), "-t", self.nodes, "-f", self.db_path, "-v", "-o", str(out_raw)]
        if info['type'] == 'paired': cmd1.extend(["-i", info['r1'], "-j", info['r2']])
        else: cmd1.extend(["-i", info['r1']])
        s1, t1, m1, _c1 = self.run_cmd(cmd1, log_f)
        if not s1: return sname, False, t1, m1, 0.0, ' '.join(cmd1)
        
        cmd2 = ["kaiju-addTaxonNames", "-t", self.nodes, "-n", self.names, "-i", str(out_raw), "-o", str(out_names)]
        s2, t2, m2, _c2 = self.run_cmd(cmd2, str(log_f)+".2")
        
        cmd3 = ["kaiju2table", "-t", self.nodes, "-n", self.names, "-r", "species", "-e", "-o", str(out_sum), str(out_raw)]
        s3, t3, m3, _c3 = self.run_cmd(cmd3, str(log_f)+".3")
        
        with open(log_f, 'a') as f1:
            for extra in [str(log_f)+".2", str(log_f)+".3"]:
                if os.path.exists(extra):
                    with open(extra, 'r') as fx: f1.write("\n\n" + fx.read())
                    Path(extra).unlink()
                    
        return sname, (s1 and s2 and s3), (t1+t2+t3), max(m1, m2, m3), 0.0, ' '.join(cmd1)

class MetabuliTool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)

        # 定义 Metabuli 特有的任务前缀，让生成的文件带有 .metabuli 标识
        job_name = f"{sname}.metabuli"
        out_rep = s_dir / f"{job_name}_report.tsv"

        # 断点续传检查对应修改为新的带前缀的文件名
        if out_rep.exists() and out_rep.stat().st_size > 0:
            return sname, True, 0.0, 0.0, 0.0, "metabuli ..."

        seq_mode = "2" if info['type'] == 'paired' else "1"
        cmd = ["metabuli", "classify"]
        if seq_mode == "2": cmd.extend([info['r1'], info['r2']])
        else: cmd.append(info['r1'])

        # 将 job_name 传给 Metabuli，它会自动生成带有该前缀的 3 个文件
        cmd.extend([self.db_path, str(s_dir), job_name, "--seq-mode", seq_mode, "--threads", str(self.threads)])

        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log")
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

class SylphTool(BaseTool):
    def process_sample(self, sname, info):
        s_dir, log_dir = self._get_sample_dirs(sname)
        out_tsv = s_dir / f"{sname}.sylph.tsv"
        if out_tsv.exists() and out_tsv.stat().st_size > 0: return sname, True, 0.0, 0.0, 0.0, "sylph ..."
        
        cmd = ["sylph", "profile", self.db_path,"--estimate-read-counts", "-t", str(self.threads),"-c","50" ]
        if info['type'] == 'paired': cmd.extend(["-1", info['r1'], "-2", info['r2']])
        else: cmd.append(info['r1'])
        
        succ, t, m, cpu_val = self.run_cmd(cmd, log_dir / f"{self.name}.log", redirect_stdout=out_tsv)
        if not succ and out_tsv.exists(): out_tsv.unlink()
        return sname, succ, t, m, cpu_val, ' '.join(cmd)

TOOL_MAP = {
    'kraken2': Kraken2Tool, 'kraken2x': Kraken2XTool, 'krakenuniq': KrakenUniqTool, 
    'centrifuger': CentrifugerTool, 'ganon': GanonTool, 'kunpeng': KunpengTool, 
    'metabuli': MetabuliTool, 'sylph': SylphTool
}

# ==========================================
# 4. 历史记录审计无损继承模块
# ==========================================
def load_existing_metrics(report_file):
    history = {}
    if Path(report_file).exists():
        try:
            with open(report_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    key = (row['Sample'], row['Tool'])
                    history[key] = row
        except Exception: pass
    return history

# ==========================================
# 5. 主流程控制：命令行解析与工业 UI 打印
# ==========================================
def main():
    try:
        import multiprocessing
        multiprocessing.set_start_method('spawn')
    except RuntimeError:
        pass

    parser = argparse.ArgumentParser(description="宏基因组统一分类器 - 完美压平版 & 工业级UI & 全能DB支持")
    
    # 核心路径参数
    parser.add_argument('-i', '--input-dir', required=True, help="FASTQ 存放目录")
    parser.add_argument('-o', '--output-dir', required=True, help="统一输出根目录 (仅存放结果)")
    parser.add_argument('--logs-dir', help="全局审计与调度日志的独立存放路径 (默认: output-dir/logs)")
    
    # 数据库配置组
    db_group = parser.add_argument_group('Database Configuration (手动或自动搜索)')
    db_group.add_argument('-d', '--db-dir', help="所有数据库统一存放的大目录 (优先手动，没指定则在此目录浅层搜索)")
    
    # 彻底暴露出所有可手动指派的工具数据库参数
    for t in list(TOOL_MAP.keys()) + ['kaiju']:
        db_group.add_argument(f'--db-{t}', help=f"【手动指定】{t.upper()} 的数据库绝对路径/前缀")
    db_group.add_argument('--kaiju-nodes', help="【手动指定】Kaiju 的 nodes.dmp 绝对路径")
    db_group.add_argument('--kaiju-names', help="【手动指定】Kaiju 的 names.dmp 绝对路径")

    # 运行配置
    parser.add_argument('--tools', nargs='+', default=['all'], help="要跑的工具集合 (默认: all)")
    parser.add_argument('--jobs', type=int, default=1, help="并行样本数 (内存消耗大时请保持 1)")
    parser.add_argument('--threads', type=int, default=30, help="单任务线程数")
    args = parser.parse_args()

    # 初始化分离式日志
    logs_dir = Path(args.logs_dir) if args.logs_dir else Path(args.output_dir) / "logs"
    logger, master_log_file = setup_logging(logs_dir)
    
    samples = find_samples(args.input_dir)
    if not samples:
        logger.error("❌ 致命错误：未找到任何 FASTQ 文件！")
        sys.exit(1)
        
    paired_cnt = sum(1 for s in samples.values() if s['type'] == 'paired')
    single_cnt = len(samples) - paired_cnt

    req_tools = list(TOOL_MAP.keys()) + ['kaiju'] if 'all' in args.tools else args.tools
    req_tools = [t.strip('{},') for t in req_tools]
    dbs = resolve_databases(args, req_tools, logger)
    
    executors = []
    skipped_tools = []
    for t in req_tools:
        if t == 'kaiju':
            if 'kaiju' in dbs and 'kaiju_nodes' in dbs:
                executors.append(KaijuTool(args.output_dir, logs_dir, args.threads, dbs))
            else: skipped_tools.append(t.upper())
        else:
            if t in dbs:
                executors.append(TOOL_MAP[t](t, args.output_dir, logs_dir, args.threads, dbs[t]))
            else: skipped_tools.append(t.upper())

    if not executors:
        logger.error("\n❌ 没有任何工具的数据库匹配成功，无法运行。")
        sys.exit(1)

    # --- 工业级美感头部打印 ---
    logger.info("="*65)
    logger.info(" 🚀 宏基因组/病毒组统一分类框架 (终极完全体) 启动 ")
    logger.info("="*65)
    logger.info(f" 📂 样本压平化输出目录: {Path(args.output_dir).absolute()}")
    logger.info(f" 📑 全局审计与报错日志: {logs_dir.absolute()}")
    logger.info(f" 📦 样本解析状态: 包含 {paired_cnt} 个双端 (PE) 样本，{single_cnt} 个单端 (SE) 样本。")
    logger.info(f" ⚡ 服务器并发配置: {args.jobs} 并发任务 × {args.threads} 线程/任务")
    logger.info("-" * 65)
    
    for ex in executors:
        logger.info(f" [就绪] {ex.name.upper().ljust(11)} | DB: {ex.db_path}")
    for st in skipped_tools:
        logger.info(f" [跳过] {st.ljust(11)} | 未手动提供且自动搜索失败。")
    logger.info("="*65)

    # 获取断点历史数据
    report_file = logs_dir / "pipeline_logs" / "sample_resource_usage.tsv"
    audit_dict = load_existing_metrics(report_file)
    
    total_tasks = len(executors) * len(samples)
    pbar = tqdm(total=total_tasks, position=0, leave=True, 
                bar_format="【全局总进度】: {l_bar}{bar}| {percentage:3.0f}% ({n_fmt}/{total_fmt}) [{postfix}]")
    global_stats = {'✅': 0, '❌': 0, '⏩': 0}
    pbar.set_postfix_str("Succ:0 | Fail:0 | Skip:0")

    task_id = 1
    for tool in executors:
        logger.info(f"\n========== 🛸 启动阶段: {tool.name.upper()} ==========")
        
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {}
            for sname, info in samples.items():
                s_type = "双端" if info['type'] == 'paired' else "单端"
                logger.info(f"[任务 {task_id:02d} | {sname}] 🚀 分发挂载 ({s_type})")
                future = pool.submit(tool.safe_process, sname, info)
                futures[future] = (task_id, sname)
                task_id += 1
                
            for future in as_completed(futures):
                tid, sname = futures[future]
                sname_ret, succ, time_sec, mem_mb, cpu_pct, cmd_str = future.result()
                
                cmd_short = cmd_str if len(cmd_str) < 100 else cmd_str[:97] + "..."
                logger.info(f"[任务 {tid:02d} | {sname}] -> {cmd_short}")
                
                key = (sname, tool.name.upper())
                
                if succ and time_sec == 0.0:
                    global_stats['⏩'] += 1
                    logger.info(f"[任务 {tid:02d} | {sname}] ⏩ 断点跳过 (结果已存在)")
                    if key not in audit_dict:
                        audit_dict[key] = {"Sample": sname, "Tool": tool.name.upper(), "Status": "Skipped", "Time(s)": "0.0", "Peak_Memory(MB)": "0.0", "CPU(%)": "0.0"}
                elif succ:
                    global_stats['✅'] += 1
                    logger.info(f"[任务 {tid:02d} | {sname}] ✅ 成功完成 (耗时: {time_sec:.2f} s | 峰值内存: {mem_mb:.2f} MB | CPU: {cpu_pct:.0f}%)")
                    audit_dict[key] = {"Sample": sname, "Tool": tool.name.upper(), "Status": "Success", "Time(s)": str(round(time_sec, 2)), "Peak_Memory(MB)": str(round(mem_mb, 2)), "CPU(%)": str(round(cpu_pct, 1))}
                else:
                    global_stats['❌'] += 1
                    logger.error(f"[任务 {tid:02d} | {sname}] ❌ 运行崩溃! 详情: {logs_dir}/{sname}/{tool.name}.log")
                    audit_dict[key] = {"Sample": sname, "Tool": tool.name.upper(), "Status": "Failed", "Time(s)": str(round(time_sec, 2)), "Peak_Memory(MB)": str(round(mem_mb, 2)), "CPU(%)": str(round(cpu_pct, 1))}
                
                pbar.update(1)
                pbar.set_postfix_str(f"Succ:{global_stats['✅']} | Fail:{global_stats['❌']} | Skip:{global_stats['⏩']}")

    pbar.close()
    
    # 无损更新 TSV 报表
    with open(report_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["Sample", "Tool", "Status", "Time(s)", "Peak_Memory(MB)", "CPU(%)"], delimiter='\t')
        writer.writeheader()
        for key in sorted(audit_dict.keys(), key=lambda x: (x[0], x[1])):
            writer.writerow(audit_dict[key])

    logger.info("\n" + "="*65)
    logger.info(" 🎉 智能批量分类管道已全部运行完毕！")
    logger.info("="*65 + "\n")

if __name__ == "__main__":
    main()
