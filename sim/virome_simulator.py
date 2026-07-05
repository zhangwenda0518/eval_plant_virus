#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
👑 Virome Benchmarker Ultimate (宏基因组基准测试 终极融合版)
特性: 
  1. 完美的 Spike-in (掺入) 策略：基于严谨的 Depth 计算靶标病毒 Reads 并掺入固定宿主背景。
  2. 长度感知丰度分配 (Length-aware Allocation)：完美还原测序仪物理原理，为 Spearman 评估提供顶级科学依据。
  3. 自动生成金标准对账单 (Manifest)：为下游 P/R/F1 和 Spearman 评估提供绝对真值。
  4. 智能任务调度 (Task-level Parallelism): 跨病毒、跨梯度大并发，彻底根除 ID 冲突 Bug。
  5. 完美断点续传 (--resume)、原子写入与 repair.sh 终极 PE 配对修复。
  6. 🛡️【终极防弹突变引擎】利用替身机制彻底绕过 mutation-simulator 的输出路径 Bug。
  7. ✨【新增】--all-in-one 模式生成 config.txt 和 SpikeIn_GroundTruth.tsv，临时文件目录规范化。
"""

import os
import sys
import time
import math
import random
import shutil
import argparse
import subprocess
import glob
import tempfile
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False

# ==========================================
# 工具 1：性能与内存监控器
# ==========================================
class PerformanceMonitor:
    def __init__(self):
        self.start_time = time.time()
        
    def report(self):
        elapsed = time.time() - self.start_time
        m, s = divmod(elapsed, 60)
        h, m = divmod(m, 60)
        time_str = f"{int(h)}小时 {int(m)}分钟 {int(s)}秒" if h > 0 else f"{int(m)}分钟 {int(s)}秒"
        
        print("\n" + "═"*50)
        print("📊 运行性能监控报告 (Performance Report)")
        print("═"*50)
        print(f"⏱️  总耗时: {time_str}")
        
        if HAS_RESOURCE:
            is_mac = sys.platform == 'darwin'
            factor = 1024 * 1024 if is_mac else 1024
            mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / factor
            mem_children = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / factor
            print(f"🧠 Python 主进程峰值内存: {mem_self:.2f} MB")
            print(f"🧬 外部工具(子进程)最高峰值内存: {mem_children:.2f} MB")
        print("═"*50 + "\n")

# ==========================================
# 工具 2：智能参数解析与文件操作
# ==========================================
def parse_number(val_str):
    s = str(val_str).strip().upper()
    try:
        if s.endswith('M'): return int(float(s[:-1]) * 1_000_000)
        elif s.endswith('K'): return int(float(s[:-1]) * 1_000)
        else: return int(float(s))
    except ValueError:
        print(f"❌ 错误: 无法将 '{val_str}' 解析为有效的数字。")
        sys.exit(1)

def parse_mut_rates(rates_list):
    result = []
    if not rates_list: return result
    for item in rates_list:
        for x in str(item).split(','):
            if x.strip(): result.append(float(x.strip()))
    return result

def resolve_targets(indir, targets_input):
    valid_exts = ('.fa', '.fasta', '.fna')
    all_files = [f for f in os.listdir(indir) if f.endswith(valid_exts)]
    if not targets_input or targets_input.lower() == 'all':
        return [os.path.join(indir, f) for f in all_files]
        
    resolved_paths = []
    for t in str(targets_input).split(','):
        t = t.strip()
        if not t: continue
        target_path = os.path.join(indir, t)
        if os.path.exists(target_path) and target_path.endswith(valid_exts):
            resolved_paths.append(target_path)
            continue
        found = False
        for ext in valid_exts:
            if os.path.exists(target_path + ext):
                resolved_paths.append(target_path + ext)
                found = True
                break
        if not found: print(f"⚠️ 警告: 找不到目标病毒 '{t}'，将跳过。")
    return resolved_paths

def is_done(files):
    if not files: return False
    return all(os.path.exists(f) and os.path.getsize(f) > 20 for f in files)

def load_config(config_path):
    config = {}
    with open(config_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                k, v = line.strip().split('\t')
                config[k] = int(v)
    return config

# ==========================================
# 工具 3：原子执行与系统命令
# ==========================================
def run_cmd(cmd):
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"\n❌ 命令执行失败:\n{' '.join(cmd)}")
        print(f"错误信息: {e.stderr.decode('utf-8')}")
        sys.exit(1)

def run_shuffle_atomic(infile, seed, outfile):
    """原子化洗牌写入：保留 .gz 后缀触发 seqkit 流式压缩"""
    if outfile.endswith(".gz"):
        tmp_out = f"{outfile}.tmp.gz"
    else:
        tmp_out = f"{outfile}.tmp"

    run_cmd(["seqkit", "shuffle", "-s", str(seed), infile, "-o", tmp_out])
    os.rename(tmp_out, outfile)

def check_dependencies(require_mut=False):
    missing = [t for t in ["art_illumina", "seqkit"] if not shutil.which(t)]
    if require_mut and not shutil.which("mutation-simulator"): missing.append("mutation-simulator")
    if missing:
        print("❌ 缺少必要的生信依赖工具: " + ", ".join(missing))
        sys.exit(1)

def concat_files_binary(input_files, output_file):
    with open(output_file, 'wb') as outfile:
        for fname in input_files:
            if not os.path.exists(fname): continue
            with open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

def get_fasta_length(fasta_path):
    return sum(len(line.strip()) for line in open(fasta_path) if not line.startswith('>'))

# ==========================================
# 核心 API 库
# ==========================================
def api_mutate(indir, outdir, rate, resume=False):
    os.makedirs(outdir, exist_ok=True)
    genomes = [f for f in os.listdir(indir) if f.endswith(('.fa', '.fasta', '.fna'))]
    path_map = {}
    print(f"  -> [Mutate] 生成突变基因组 (Rate: {rate}%)...")
    
    skipped = 0
    for genome in genomes:
        in_path = os.path.join(indir, genome)
        base = os.path.splitext(genome)[0]
        final_path = os.path.join(outdir, f"{base}_mut_{rate}pct.fa")
        
        path_map[in_path] = final_path
        
        if resume and os.path.exists(final_path) and os.path.getsize(final_path) > 10:
            # 验证是有效 FASTA（排除 .fai / .txt 等垃圾文件）
            with open(final_path) as check_f:
                first_line = check_f.readline().strip()
            if first_line.startswith('>'):
                skipped += 1
                continue
            else:
                os.remove(final_path)  # 无效文件，删除并重新生成
            
        tmp_in = os.path.join(outdir, f"tmp_in_{genome}")
        shutil.copy(in_path, tmp_in)
        
        # 突变模拟
        cmd = ["mutation-simulator", tmp_in, "args", "-sn", str(rate/100.0), "-in", "0.001"]
        run_cmd(cmd)
        
        tmp_base = os.path.splitext(f"tmp_in_{genome}")[0]
        expected_out = os.path.join(outdir, f"{tmp_base}_mutated.fasta")
        
        # 防弹机制：抓取任何可能的输出名称
        if os.path.exists(expected_out):
            shutil.move(expected_out, final_path)
        else:
            found = glob.glob(os.path.join(outdir, f"{tmp_base}*"))
            fasta_only = [f for f in found if f.endswith(('.fa','.fasta')) and 'mutated' in os.path.basename(f)]
            if fasta_only:
                shutil.move(fasta_only[0], final_path)
            elif found:
                for f in sorted(found):
                    try:
                        with open(f) as cf:
                            if cf.readline().strip().startswith('>'):
                                shutil.move(f, final_path)
                                print(f"  ⚠️ 防弹恢复: {os.path.basename(f)} -> {os.path.basename(final_path)}")
                                break
                    except: pass
                else:
                    print(f"  ⚠️ mutation-simulator 对 {genome} 静默失败，跳过")
            else:
                print(f"  ⚠️ mutation-simulator 对 {genome} 静默失败，跳过")
                
        if os.path.exists(tmp_in): os.remove(tmp_in)
        txt_log = os.path.join(outdir, f"{tmp_base}_mutated_mutations.txt")
        if os.path.exists(txt_log): os.remove(txt_log)
            
    if skipped > 0: print(f"  ⏭️ [断点续传] 已跳过 {skipped} 个完整存在的突变基因组")
    return path_map

def api_gen_config(indir, read_len, depth=None, total_reads=None):
    """👑 核心修改: 引入长度感知的摩尔丰度分配 (Length-aware Molar Allocation)"""
    genomes = [os.path.join(indir, f) for f in os.listdir(indir) if f.endswith(('.fa', '.fasta', '.fna'))]
    config_dict = {}
    
    if depth:
        # 模式1: 统一组装深度模式
        for g in genomes:
            r = math.ceil((depth * get_fasta_length(g)) / read_len)
            config_dict[g] = r + 1 if r % 2 != 0 else r
            
    elif total_reads:
        # 模式2: Spearman丰度测序模式 (符合真实测序仪物理原理)
        molar_abundances = [random.lognormvariate(0, 1.5) for _ in range(len(genomes))]
        lengths = [get_fasta_length(g) for g in genomes]
        
        # 测序 Reads 数正比于 (摩尔丰度 * 基因组长度)
        seq_weights = [molar_abundances[i] * lengths[i] for i in range(len(genomes))]
        total_weight = sum(seq_weights)
        norm_weights = [w / total_weight for w in seq_weights]
        
        allocated = 0
        for i, g in enumerate(genomes):
            if i == len(genomes) - 1:
                r = total_reads - allocated
            else:
                r = int(total_reads * norm_weights[i])
                
            r = r + 1 if r % 2 != 0 else r
            config_dict[g] = r
            allocated += r
            
    return config_dict

def api_run_sim(config_dict, out_prefix, mode, read_len, profile, seed, threads=1, frag_mean=250, frag_std=15, resume=False):
    final_r1 = f"{out_prefix}_PE_R1.fastq.gz" if mode == "PE" else f"{out_prefix}_SE.fastq.gz"
    final_r2 = f"{out_prefix}_PE_R2.fastq.gz" if mode == "PE" else None
    
    if resume and is_done([final_r1, final_r2] if mode == "PE" else [final_r1]):
        print(f"  ⏭️ [断点续传] 数据已存在，跳过测序: {os.path.basename(final_r1)}")
        return final_r1, final_r2

    tmp_dir = f"{out_prefix}_tmp_ART"
    os.makedirs(tmp_dir, exist_ok=True)
    
    def run_art_genome(task_args):
        ref, fold_cov, idx = task_args
        base = os.path.splitext(os.path.basename(ref))[0]
        out_base = os.path.join(tmp_dir, f"{base}_sim{idx}_")
        
        art_seed = int(seed + idx * 777) % 2147483647 
        
        cmd = ["art_illumina", "-ss", profile, "-i", ref, "-l", str(read_len), "-na", "-q", "-qs", "40", "-qs2", "40", "-rs", str(art_seed), "-o", out_base]
        if mode == "PE": 
            cmd.extend(["-p", "-m", str(frag_mean), "-s", str(frag_std), "-f", f"{fold_cov:.6f}"])
        else: 
            cmd.extend(["-f", f"{fold_cov:.6f}"])
        run_cmd(cmd)

    tasks = []
    task_idx = 0
    for ref, target_reads in config_dict.items():
        if target_reads <= 0: continue
        genome_len = get_fasta_length(ref)
        if genome_len == 0: continue
        total_fold_cov = (target_reads * read_len) / float(genome_len)
        tasks.append((ref, total_fold_cov, task_idx))
        task_idx += 1

    active_threads = min(threads, len(tasks) if tasks else 1)
    with ThreadPoolExecutor(max_workers=active_threads) as executor:
        list(executor.map(run_art_genome, tasks))

    # 合并与洗牌
    if mode == "PE":
        raw_r1, raw_r2 = os.path.join(tmp_dir, "raw_R1.fq"), os.path.join(tmp_dir, "raw_R2.fq")
        concat_files_binary(glob.glob(os.path.join(tmp_dir, "*1.fq")), raw_r1)
        concat_files_binary(glob.glob(os.path.join(tmp_dir, "*2.fq")), raw_r2)
        tmp_r1, tmp_r2 = f"{out_prefix}_tmpR1.fq.gz", f"{out_prefix}_tmpR2.fq.gz"
        run_shuffle_atomic(raw_r1, seed, tmp_r1)
        run_shuffle_atomic(raw_r2, seed, tmp_r2)
        
        if shutil.which("repair.sh"):
            run_cmd(["repair.sh", f"in1={tmp_r1}", f"in2={tmp_r2}", f"out1={final_r1}", f"out2={final_r2}", "overwrite=true", "qin=33"])
            os.remove(tmp_r1); os.remove(tmp_r2)
        else:
            os.rename(tmp_r1, final_r1); os.rename(tmp_r2, final_r2)
    else:
        raw_se = os.path.join(tmp_dir, "raw_SE.fq")
        concat_files_binary(glob.glob(os.path.join(tmp_dir, "*.fq")), raw_se)
        run_shuffle_atomic(raw_se, seed, final_r1)
        
    shutil.rmtree(tmp_dir)
    return final_r1, final_r2

def api_subsample(r1, r2, se, mode, outdir, depths, repeats, threads, resume=False):
    os.makedirs(outdir, exist_ok=True)
    print(f"  -> [抽样] 启动 Jackknife 多深度并发压缩抽样...")
    
    if mode == "PE":
        base_name = os.path.basename(r1).replace("_PE_R1.fastq.gz", "").replace("_R1.fastq.gz", "")
    else:
        base_name = os.path.basename(se).replace("_SE.fastq.gz", "").replace(".fastq.gz", "")
        
    def sample_single(infile, depth, seed, outfile):
        if resume and is_done([outfile]): return 1
        if outfile.endswith(".gz"): tmp_out = f"{outfile}.tmp.gz"
        else: tmp_out = f"{outfile}.tmp"
            
        run_cmd(["seqkit", "sample", "-s", str(seed), "-n", str(depth), infile, "-o", tmp_out])
        os.rename(tmp_out, outfile)
        return 0

    futures = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for d in depths:
            for s in range(1, repeats + 1):
                if mode == "PE":
                    futures.append(executor.submit(sample_single, r1, d, s, os.path.join(outdir, f"{base_name}_sub_{d}_rep{s}_R1.fastq.gz")))
                    futures.append(executor.submit(sample_single, r2, d, s, os.path.join(outdir, f"{base_name}_sub_{d}_rep{s}_R2.fastq.gz")))
                else:
                    futures.append(executor.submit(sample_single, se, d, s, os.path.join(outdir, f"{base_name}_sub_{d}_rep{s}_SE.fastq.gz")))
                    
    skipped = sum(f.result() for f in futures)
    if skipped > 0: print(f"  ⏭️ [断点续传] 跳过了 {skipped} 个已完成的抽样文件")

def api_lod_test(bg_ref, target_paths, bg_reads, depths, read_len, mode, outdir, threads, seed, resume=False, all_in_one=False):
    """
    💎 最核心的 Spike-in (掺入) 策略：固定宿主背景，按指定的 Depth 掺入病毒 Reads，并生成金标准报告。
    改进: all-in-one 模式输出 config.txt 和 SpikeIn_GroundTruth.tsv，所有临时文件放入临时目录。
    """
    os.makedirs(outdir, exist_ok=True)
    
    # 创建临时根目录，存放所有中间文件
    temp_root = tempfile.mkdtemp(prefix="lod_temp_", dir=outdir)
    try:
        print(f"\n[LoD 测试] 正在生成共享宿主背景 (指定 {bg_reads} Reads，单线程运行约需数分钟)...")
        bg_reads_even = bg_reads + 1 if bg_reads % 2 != 0 else bg_reads
        
        bg_prefix = os.path.join(temp_root, "Shared_Background")
        tmp_bg_r1, tmp_bg_r2 = api_run_sim({bg_ref: bg_reads_even}, bg_prefix, mode, read_len, "HS25", seed, threads=2, resume=resume)
        
        if all_in_one:
            # ---------- all-in-one 模式：混合所有病毒到一个样本，并生成 config.txt 和真值表 ----------
            out_prefix = os.path.join(outdir, "Host_Depletion_Mixed")
            out_r1 = f"{out_prefix}_PE_R1.fastq.gz"
            out_r2 = f"{out_prefix}_PE_R2.fastq.gz"
            if resume and is_done([out_r1, out_r2]):
                print("  [resume] All-in-one exists, skip")
                return
            
            if len(depths) != 1:
                print("⚠️ Warning: --all-in-one 模式建议只给一个 depth 值，将使用 depths[0]")
            target_depth = depths[0]
            
            # 收集每个病毒的配置信息，同时写入 config.txt 和真值表
            config_lines = []  # 用于 config.txt: path\treads
            truth_lines = []   # 用于 SpikeIn_GroundTruth.tsv: 病毒名、深度、长度、Host_Reads、Virus_Reads、Total_Reads、丰度
            truth_lines.append("Virus_File\tTarget_Depth(x)\tGenome_Length(bp)\tHost_Reads\tVirus_Reads\tTotal_Reads\tTrue_Virus_Abundance(%)")
            
            all_virus_r1, all_virus_r2 = [], []
            for v_path in tqdm(target_paths, desc="Generating virus reads"):
                gl = get_fasta_length(v_path)
                if gl == 0:
                    continue
                raw_v_reads = (target_depth * gl) / read_len
                v_reads = math.ceil(raw_v_reads)
                v_reads = v_reads + 1 if v_reads % 2 != 0 else v_reads
                if v_reads <= 0:
                    continue
                
                # config.txt 格式：路径\treads
                config_lines.append(f"{v_path}\t{v_reads}")
                
                total = bg_reads_even + v_reads
                abundance = (v_reads / total) * 100
                truth_lines.append(f"{os.path.basename(v_path)}\t{target_depth}\t{gl}\t{bg_reads_even}\t{v_reads}\t{total}\t{abundance:.6f}")
                
                v_prefix = os.path.join(temp_root, f"virus_{os.path.basename(v_path)}")
                v1, v2 = api_run_sim({v_path: v_reads}, v_prefix, mode, read_len, "HS25", seed, threads=2, resume=resume)
                all_virus_r1.append(v1)
                if v2:
                    all_virus_r2.append(v2)
            
            # 写入 config.txt
            config_path = os.path.join(outdir, "config.txt")
            with open(config_path, 'w') as f:
                f.write("\n".join(config_lines))
            print(f"✅ 病毒掺入配置已保存（config.txt 格式）: {config_path}")
            
            # 写入真值表 SpikeIn_GroundTruth.tsv
            truth_path = os.path.join(outdir, "SpikeIn_GroundTruth.tsv")
            with open(truth_path, 'w') as f:
                f.write("\n".join(truth_lines))
            print(f"✅ 病毒掺入真值表已保存: {truth_path}")
            
            # 合并、洗牌、修复
            raw_r1 = os.path.join(temp_root, "mixed_raw_R1.fq.gz")
            raw_r2 = os.path.join(temp_root, "mixed_raw_R2.fq.gz")
            concat_files_binary([tmp_bg_r1] + all_virus_r1, raw_r1)
            concat_files_binary([tmp_bg_r2] + all_virus_r2, raw_r2)
            
            tmp_r1 = os.path.join(temp_root, "mixed_tmp_R1.fq.gz")
            tmp_r2 = os.path.join(temp_root, "mixed_tmp_R2.fq.gz")
            run_shuffle_atomic(raw_r1, seed, tmp_r1)
            run_shuffle_atomic(raw_r2, seed, tmp_r2)
            
            if shutil.which("repair.sh"):
                run_cmd(["repair.sh", f"in1={tmp_r1}", f"in2={tmp_r2}", f"out1={out_r1}", f"out2={out_r2}", "overwrite=true", "qin=33"])
            else:
                os.rename(tmp_r1, out_r1)
                os.rename(tmp_r2, out_r2)
            
            print(f"✅ All-in-one 混合样本生成: {out_r1} 和 {out_r2}")
            return  # all-in-one 结束，不执行下面的逐任务循环
        
        # ---------- 非 all-in-one 模式：每个病毒/深度独立样本（原逻辑，仅临时目录调整） ----------
        tasks = [(v_path, d) for v_path in target_paths for d in depths]
        print(f"  -> 准备完毕！启动 {threads} 个并发进程，跨 {len(target_paths)} 个病毒与 {len(depths)} 个深度梯度同时掺入...")
        
        def process_lod_task(task_args):
            v_path, target_depth = task_args
            virus_name = os.path.splitext(os.path.basename(v_path))[0]
            mix_prefix = os.path.join(outdir, f"LoD_Mixed_{virus_name}_{target_depth}x")
            mix_r1 = f"{mix_prefix}_PE_R1.fastq.gz" if mode == "PE" else f"{mix_prefix}_SE.fastq.gz"
            mix_r2 = f"{mix_prefix}_PE_R2.fastq.gz" if mode == "PE" else None
            
            if resume and is_done([mix_r1, mix_r2] if mode == "PE" else [mix_r1]):
                return
            
            genome_len = get_fasta_length(v_path)
            if genome_len == 0:
                return
            
            v_reads = math.ceil((target_depth * genome_len) / read_len)
            v_reads = v_reads + 1 if v_reads % 2 != 0 else v_reads
            if v_reads <= 0:
                return
            
            # 病毒模拟放入临时目录
            v_prefix = os.path.join(temp_root, f"tmp_{virus_name}_{target_depth}x")
            art_seed = int(seed + target_depth * 10000) % 2147483647
            tmp_v_r1, tmp_v_r2 = api_run_sim({v_path: v_reads}, v_prefix, mode, read_len, "HS25", art_seed, threads=2, resume=resume)
            
            if mode == "PE":
                raw_r1 = os.path.join(temp_root, f"{virus_name}_{target_depth}x_raw_R1.fq.gz")
                raw_r2 = os.path.join(temp_root, f"{virus_name}_{target_depth}x_raw_R2.fq.gz")
                concat_files_binary([tmp_bg_r1, tmp_v_r1], raw_r1)
                concat_files_binary([tmp_bg_r2, tmp_v_r2], raw_r2)
                
                tmp_r1 = os.path.join(temp_root, f"{virus_name}_{target_depth}x_tmp_R1.fq.gz")
                tmp_r2 = os.path.join(temp_root, f"{virus_name}_{target_depth}x_tmp_R2.fq.gz")
                run_shuffle_atomic(raw_r1, seed, tmp_r1)
                run_shuffle_atomic(raw_r2, seed, tmp_r2)
                
                if shutil.which("repair.sh"):
                    run_cmd(["repair.sh", f"in1={tmp_r1}", f"in2={tmp_r2}", f"out1={mix_r1}", f"out2={mix_r2}", "overwrite=true", "qin=33"])
                else:
                    os.rename(tmp_r1, mix_r1)
                    os.rename(tmp_r2, mix_r2)
                # 清理病毒临时文件
                for f in [raw_r1, raw_r2, tmp_v_r1, tmp_v_r2]:
                    if os.path.exists(f):
                        os.remove(f)
            else:
                # SE 模式类似处理
                raw_se = os.path.join(temp_root, f"{virus_name}_{target_depth}x_raw_SE.fq.gz")
                concat_files_binary([tmp_bg_r1, tmp_v_r1], raw_se)
                run_shuffle_atomic(raw_se, seed, mix_r1)
                os.remove(raw_se)
                os.remove(tmp_v_r1)
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            list(tqdm(executor.map(process_lod_task, tasks), total=len(tasks), desc="并行混样进度"))
        
        # 生成金标准对账单（原逻辑）
        manifest_path = os.path.join(outdir, "LoD_GroundTruth_Manifest.tsv")
        print(f"\n  -> [收尾] 正在生成金标准对账单: {manifest_path}")
        with open(manifest_path, 'w') as f:
            f.write("Sample_Name\tTarget_Virus\tTarget_Depth(x)\tHost_Reads\tVirus_Reads\tTotal_Reads\tTrue_Virus_Abundance(%)\n")
            for v_path in target_paths:
                genome_len = get_fasta_length(v_path)
                if genome_len == 0:
                    continue
                virus_name = os.path.splitext(os.path.basename(v_path))[0]
                for d in depths:
                    v_reads = math.ceil((d * genome_len) / read_len)
                    v_reads = v_reads + 1 if v_reads % 2 != 0 else v_reads
                    if v_reads <= 0:
                        continue
                    total = bg_reads_even + v_reads
                    true_abundance = (v_reads / total) * 100 if total > 0 else 0
                    sample_name = f"LoD_Mixed_{virus_name}_{d}x"
                    f.write(f"{sample_name}\t{virus_name}\t{d}\t{bg_reads_even}\t{v_reads}\t{total}\t{true_abundance:.6f}\n")
        print(f"✅ 金标准对账单已保存，下游可直接用于 P/R/F1 和 Spearman 评估！")
    
    finally:
        # 清理所有临时文件
        shutil.rmtree(temp_root, ignore_errors=True)
        print(f"🧹 已清理临时目录: {temp_root}")

# ==========================================
# 独立子命令 Wrapper
# ==========================================
def cmd_mutate(args):
    check_dependencies(require_mut=True)
    for r in parse_mut_rates(args.rate): api_mutate(args.indir, args.outdir, r, args.resume)
    print(f"✅ 突变生成完毕: {args.outdir}")

def cmd_gen_config(args):
    print("🌟 开始生成群落丰度配置文件...")
    if args.resume and os.path.exists(args.outconfig):
        print(f"⏭️ [断点续传] 配置文件已存在: {args.outconfig}")
        return
    config = api_gen_config(args.indir, args.read_len, parse_number(args.depth) if args.depth else None, parse_number(args.total_reads) if args.total_reads else None)
    with open(args.outconfig, 'w') as f:
        for g, r in config.items(): f.write(f"{g}\t{r}\n")
    print(f"✅ 配置文件已保存: {args.outconfig}")

def cmd_simulator(args):
    check_dependencies()
    print("🚀 开始多线程模拟测序并洗牌...")
    api_run_sim(load_config(args.config), args.out, args.mode, args.read_len, args.profile, args.seed, args.threads, args.frag_mean, args.frag_std, args.resume)
    print("✅ 模拟与大洗牌完成！")

def cmd_subsample(args):
    check_dependencies()
    api_subsample(args.r1, args.r2, args.se, args.mode, args.outdir, [parse_number(d) for d in args.depths], args.repeats, args.threads, args.resume)
    print(f"✅ 子采样全部完成: {args.outdir}")

def cmd_lod_mix(args):
    check_dependencies()
    print("🦠 启动独立 LoD 大海捞针测试...")
    target_paths = resolve_targets(args.indir, args.targets)
    api_lod_test(args.bgref, target_paths, args.bg_reads, args.depths, args.read_len, args.mode, args.outdir, args.threads, args.seed, args.resume, args.all_in_one)
    print("✅ LoD 测试数据集生成完毕！")

def cmd_benchmark(args):
    check_dependencies(require_mut=True)
    os.makedirs(args.outdir, exist_ok=True)
    print("\n" + "="*50)
    print("🚀 启动 Ultimate Benchmark 终极流水线")
    print("="*50)

    base_config_path = os.path.join(args.outdir, "Base_Abundance_Profile.txt")
    if args.resume and os.path.exists(base_config_path):
        print("\n[1/3] ⏭️ [断点续传] 读取已存的全局黄金丰度比例...")
        base_config = load_config(base_config_path)
    else:
        print("\n[1/3] 🎲 锁定全局标准群落丰度 (变量控制)...")
        if args.depth:
            base_config = api_gen_config(args.indir, args.read_len, depth=parse_number(args.depth))
        else:
            base_config = api_gen_config(args.indir, args.read_len, total_reads=args.total_reads)
        with open(base_config_path, 'w') as f:
            for g, r in base_config.items(): f.write(f"{g}\t{r}\n")

    rates = parse_mut_rates(args.mut_rates)
    depths = [parse_number(d) for d in args.depths]
    
    for rate in rates:
        rate_dir = os.path.join(args.outdir, f"Dataset_Mut_{rate}pct")
        os.makedirs(rate_dir, exist_ok=True)
        
        current_config = base_config.copy()
        if rate != 0:
            path_map = api_mutate(args.indir, os.path.join(rate_dir, "Mutated_Genomes"), rate, args.resume)
            current_config = {path_map[orig]: reads for orig, reads in base_config.items()}
                
        m_prefix = os.path.join(rate_dir, f"Master_{rate}pct")
        r1, r2 = api_run_sim(current_config, m_prefix, args.mode, args.read_len, args.profile, args.seed, args.threads, resume=args.resume)
        
        sub_dir = os.path.join(rate_dir, "Jackknife_Subsamples")
        if args.mode == "PE": api_subsample(r1, r2, None, "PE", sub_dir, depths, args.repeats, args.threads, args.resume)
        else: api_subsample(None, None, r1, "SE", sub_dir, depths, args.repeats, args.threads, args.resume)

    if args.bgref and os.path.exists(args.bgref):
        lod_dir = os.path.join(args.outdir, "Dataset_LoD_Test")
        target_paths = resolve_targets(args.indir, args.targets)
        # 宿主背景 10M reads
        api_lod_test(args.bgref, target_paths, 10000000, [0.1, 0.5, 1.0, 5.0, 10.0, 50.0], args.read_len, args.mode, lod_dir, args.threads, args.seed, args.resume)
    
    print("\n" + "="*50)
    print(f"🎉 终极流水线执行完毕！结果保存在: {os.path.abspath(args.outdir)}")
    print("="*50)

# ==========================================
# CLI 解析引擎
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="👑 Virome Simulator Ultimate (全功能宏基因组数据生成引擎)")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common_args(p):
        p.add_argument("--resume", action="store_true", help="启用断点续传")

    # Mutate Command
    p_mut = subparsers.add_parser("mutate")
    p_mut.add_argument("-i", "--indir", required=True); p_mut.add_argument("-o", "--outdir", required=True)
    p_mut.add_argument("-r", "--rate", nargs='+', default=["5"]); add_common_args(p_mut); p_mut.set_defaults(func=cmd_mutate)

    # Gen Config Command
    p_cfg = subparsers.add_parser("gen-config")
    p_cfg.add_argument("-i", "--indir", required=True); p_cfg.add_argument("-o", "--outconfig", default="config.txt")
    p_cfg.add_argument("-l", "--read-len", type=int, default=150)
    g = p_cfg.add_mutually_exclusive_group(required=True)
    g.add_argument("-d", "--depth", type=str); g.add_argument("-t", "--total-reads", type=str)
    add_common_args(p_cfg); p_cfg.set_defaults(func=cmd_gen_config)

    # Simulator Command
    p_sim = subparsers.add_parser("simulator")
    p_sim.add_argument("-c", "--config", required=True); p_sim.add_argument("-o", "--out", default="Simulated_Data")
    p_sim.add_argument("--mode", choices=["SE", "PE"], default="PE"); p_sim.add_argument("--threads", type=int, default=8)
    p_sim.add_argument("-l", "--read-len", type=int, default=150); p_sim.add_argument("--profile", default="HS25")
    p_sim.add_argument("--seed", type=int, default=42); 
    p_sim.add_argument("-s", "--frag-std", type=int, default=50); add_common_args(p_sim); p_sim.set_defaults(func=cmd_simulator)
    p_sim.add_argument("-m", "--frag-mean", type=int, default=300)

    # Subsample Command
    p_sub = subparsers.add_parser("subsample")
    p_sub.add_argument("--mode", choices=["SE", "PE"], required=True); p_sub.add_argument("--r1", help="R1")
    p_sub.add_argument("--r2", help="R2"); p_sub.add_argument("--se", help="SE"); p_sub.add_argument("-o", "--outdir", default="Subsampled_Data")
    p_sub.add_argument("-d", "--depths", nargs='+', required=True); p_sub.add_argument("-r", "--repeats", type=int, default=5)
    p_sub.add_argument("--threads", type=int, default=8); add_common_args(p_sub); p_sub.set_defaults(func=cmd_subsample)

    # LoD Mix Command (Depth-based Spike-in)
    p_lod = subparsers.add_parser("LoD_mix")
    p_lod.add_argument("-i", "--indir", required=True); p_lod.add_argument("--bgref", required=True)
    p_lod.add_argument("--targets", default="all"); p_lod.add_argument("--bg-reads", type=parse_number, default=1000000)
    p_lod.add_argument("--depths", nargs='+', type=float, default=[0.5,1.0,5.0, 10.0,20.0, 50.0,100.0,200.0,500.0,1000.0])
    p_lod.add_argument("--mode", choices=["SE", "PE"], default="PE")
    p_lod.add_argument("-l", "--read-len", type=int, default=150); p_lod.add_argument("-o", "--outdir", default="LoD_Dataset")
    p_lod.add_argument("--seed", type=int, default=42); p_lod.add_argument("--threads", type=int, default=8)
    p_lod.add_argument("--all-in-one", action="store_true", help="合并所有病毒到1个样本，并生成config.txt和SpikeIn_GroundTruth.tsv")
    add_common_args(p_lod); p_lod.set_defaults(func=cmd_lod_mix)

    # Benchmark Command
    p_bench = subparsers.add_parser("benchmark")
    p_bench.add_argument("-i", "--indir", required=True); p_bench.add_argument("-o", "--outdir", default="Ultimate_Benchmark_Results")
    p_bench.add_argument("-t", "--total-reads", type=parse_number, default=20000000); p_bench.add_argument("--mode", choices=["SE", "PE"], default="PE")
    p_bench.add_argument("--threads", type=int, default=8); p_bench.add_argument("--mut-rates", nargs='+', default=["0", "5", "15"])
    p_bench.add_argument("--depths", nargs='+', default=["50k", "250k", "1M"]); p_bench.add_argument("--repeats", type=int, default=5)
    p_bench.add_argument("--bgref", default=""); p_bench.add_argument("--targets", default="all")
    p_bench.add_argument("--depth", type=str, default=None, help="统一覆盖深度(如200)，均匀分配reads（替代--total-reads随机分配）")
    p_bench.add_argument("-l", "--read-len", type=int, default=150); p_bench.add_argument("--profile", default="HS25")
    p_bench.add_argument("--seed", type=int, default=42); add_common_args(p_bench); p_bench.set_defaults(func=cmd_benchmark)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    monitor = PerformanceMonitor()
    try:
        main()
    finally:
        monitor.report()
