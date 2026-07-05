#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Virus Assembly Benchmark — Phase 2: 嵌合 Contig 检测
================================================================================

功能:
  1. 为每个 Group 建 BLAST 数据库
  2. 并行 BLASTn + detect_chimeric_contigs.py
  3. 解析嵌合报告 → chimeric_summary.tsv

断点续跑: 检查 benchmark_results/chimeric/{sample}_{tool}_chimeric.tsv 是否存在

用法:
  python scripts/benchmark_chimeric.py -i step6_assemblies -r step6_ref_viruses.fasta -o benchmark_results --mode 7 -j 8 -t 4

输出:
  - benchmark_results/blast_db/group_N/            ← 每个 Group 的 BLAST 数据库
  - benchmark_results/chimeric/{sample}_{tool}_*.tsv
  - benchmark_results/chimeric_summary.tsv
"""

import os, sys, argparse, subprocess, shutil, gzip
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import (
    discover_samples, get_assembly_files, parse_dir_name, get_existing_file,
    FILE_PATTERNS
)
from benchmark_mq import _split_fasta_by_group


def count_contigs_in_fasta(fasta_path):
    """统计 FASTA 中的 contig 数量"""
    count = 0
    try:
        handler = gzip.open if fasta_path.endswith('.gz') else open
        with handler(fasta_path, 'rt' if fasta_path.endswith('.gz') else 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    except Exception:
        return 0


def count_chimeric(report_path):
    """统计嵌合报告中的嵌合 contig 数"""
    if not os.path.exists(report_path):
        return 0
    try:
        with open(report_path, 'r') as f:
            return max(0, len(f.readlines()) - 1)
    except Exception:
        return 0


def run_blastn_and_detect(args_tuple):
    """单任务: BLASTn + detect_chimeric"""
    ctg_path, tool, sn, out_tsv, rpt_tsv, db_name, chimeric_script, threads = args_tuple

    # 1. BLASTn
    blast_cmd = [
        "blastn", "-query", ctg_path, "-db", db_name,
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-out", out_tsv, "-evalue", "1e-5", "-num_threads", str(threads),
        "-max_target_seqs", "10",
    ]
    try:
        subprocess.run(blast_cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=7200)
    except Exception as e:
        return f"❌ BLASTn {sn}/{tool}: {e}"

    # 2. detect_chimeric: -i 输入 --out 输出 (注意: -o 是 max_overlap)
    detect_cmd = [sys.executable, chimeric_script, "-i", out_tsv, "--out", rpt_tsv]
    try:
        subprocess.run(detect_cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=600)
    except Exception as e:
        return f"❌ Chimeric {sn}/{tool}: {e}"
    return f"✅ {sn}/{tool}"


def main():
    parser = argparse.ArgumentParser(description="🧬 Phase 2: 嵌合 Contig 检测")
    parser.add_argument("-i", "--input-dir", default="step6_assemblies")
    parser.add_argument("-r", "--ref", required=True,
                        help="参考病毒: 合并 FASTA 或包含病毒序列的目录")
    parser.add_argument("-o", "--out-dir", default="benchmark_results")
    parser.add_argument("--chimeric-script", default=None,
                        help="detect_chimeric_contigs.py 脚本路径 (默认: 同目录)")
    parser.add_argument("--mode", choices=['4', '6', '7', '10'], default='6')
    parser.add_argument("-j", "--jobs", type=int, default=8, help="并行 BLASTn 任务数")
    parser.add_argument("-t", "--threads", type=int, default=4, help="单个 BLASTn 线程数")
    parser.add_argument("--force", action="store_true", help="强制重新运行")
    
    # 🌟 新增：智能控制参考序列是否按照 group 拆分的参数
    parser.add_argument("--split-ref", choices=['auto', 'always', 'never'], default='auto',
                        help="处理合并的 FASTA 时是否按 Group 拆分 (默认 auto: 根据样本自动判断)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ---- 0. 检查依赖 ----
    if args.chimeric_script is None:
        args.chimeric_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "detect_chimeric_contigs.py")
    blastn_exe = shutil.which("blastn")
    makeblastdb_exe = shutil.which("makeblastdb")
    if not blastn_exe or not makeblastdb_exe:
        print("❌ 未找到 BLAST+！请安装: conda install -c bioconda blast")
        return 1
    if not os.path.exists(args.chimeric_script):
        print(f"❌ 未找到: {args.chimeric_script}")
        return 1
    print(f"✅ blastn: {blastn_exe}")

    # ---- 1. 提前扫描样本 (用于自动判断是否需要拆分) ----
    samples = discover_samples(args.input_dir)
    if not samples:
        print("❌ 未发现任何样本目录！")
        return 1
    
    # 判断样本中是否具有 group 信息
    has_groups = any(s.get('group', 0) != 0 for s in samples)

    # ---- 1.5 准备 Group 参考序列 ----
    ref_path = os.path.abspath(args.ref)
    group_ref_map = {}

    if os.path.isfile(ref_path):
        # 决定拆分策略
        do_group_split = False
        if args.split_ref == 'always':
            do_group_split = True
        elif args.split_ref == 'auto' and has_groups:
            do_group_split = True

        if do_group_split:
            print(f"📦 输入为合并 FASTA: {ref_path}，[策略: {args.split_ref}] 正在按 Group 拆分...")
            ref_dir_base = os.path.join(args.out_dir, "_ref_split_grouped")
            group_ref_map = _split_fasta_by_group(ref_path, ref_dir_base)
            print(f"📁 拆分完成: {len(group_ref_map)} 个 group 参考目录")
        else:
            print(f"📦 输入为合并 FASTA: {ref_path}，[策略: {args.split_ref}] 不分组，但正在按病毒拆分为独立文件...")
            flat_ref_dir = os.path.join(args.out_dir, "_ref_split_flat")
            os.makedirs(flat_ref_dir, exist_ok=True)
            
            # 高效流式读取，将合并序列拆成一个个独立文件
            with open(ref_path, 'r') as fh:
                current_id, current_seq = None, []
                for line in fh:
                    if line.startswith('>'):
                        if current_id:
                            with open(os.path.join(flat_ref_dir, f"{current_id}.fasta"), 'w') as out:
                                out.write(f">{current_id}\n{''.join(current_seq)}")
                        current_id = line.lstrip('>').strip().split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_id:
                    with open(os.path.join(flat_ref_dir, f"{current_id}.fasta"), 'w') as out:
                        out.write(f">{current_id}\n{''.join(current_seq)}")
            
            # 统一分派为 group 0 (全局比对)
            group_ref_map = {0: flat_ref_dir}
            print(f"📁 独立病毒拆分完成，已全部放入全局参考目录。")

    elif os.path.isdir(ref_path):
        for gdir in sorted(os.listdir(ref_path)):
            gpath = os.path.join(ref_path, gdir)
            if os.path.isdir(gpath) and gdir.startswith("group_"):
                try:
                    gnum = int(gdir.split("_")[1])
                    group_ref_map[gnum] = gpath
                except ValueError:
                    continue
        
        # 兼容性修复：如果目录里全是扁平文件，没有 group_ 文件夹，把整个目录当成全局 Group 0
        if not group_ref_map:
            print(f"⚠️ 参考目录中未找到 group_ 子目录，将整个目录视为全局参考 (Group 0)")
            group_ref_map = {0: ref_path}
        else:
            print(f"📁 参考目录: {len(group_ref_map)} 个 group")
    else:
        print(f"❌ 参考路径不存在: {ref_path}")
        return 1

    # ---- 2. 为每个 Group 建 BLAST 数据库 ----
    blast_db_base = os.path.join(args.out_dir, "blast_db")
    os.makedirs(blast_db_base, exist_ok=True)
    group_db_map = {}

    for g, gdir in group_ref_map.items():
        # 合并该 group 的 FASTA
        merged = os.path.join(blast_db_base, f"group_{g}.fasta")
        if not os.path.exists(merged):
            with open(merged, 'w') as mf:
                for fname in sorted(os.listdir(gdir)):
                    if fname.endswith('.fasta'):
                        with open(os.path.join(gdir, fname), 'r') as inf:
                            mf.write(inf.read())

        db_name = os.path.join(blast_db_base, f"group_{g}")
        if not os.path.exists(db_name + ".nin"):
            print(f"📦 建 BLAST 数据库 group_{g} ({len([f for f in os.listdir(gdir) if f.endswith('.fasta')])} viruses)")
            subprocess.run(["makeblastdb", "-in", merged, "-dbtype", "nucl",
                            "-out", db_name, "-parse_seqids"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        group_db_map[g] = db_name
    print(f"📦 {len(group_db_map)} 个 BLAST 数据库就绪")

    # ---- 3. 收集 BLASTn 任务 ----
    chim_dir = os.path.join(args.out_dir, "chimeric")
    os.makedirs(chim_dir, exist_ok=True)
    blast_tasks = []
    skipped_no_ref = 0
    skipped_done = 0

    for s in samples:
        g = s.get('group', 0)
        db_name = group_db_map.get(g)
        if db_name is None:
            skipped_no_ref += 1
            continue
        files = get_assembly_files(s['sample_dir'], s['base_name'], args.mode)
        for tool, fpath in files:
            sn = s['base_name']
            # 加上 group 前缀，避免不同 group 同名 sample 互相覆盖
            prefix = f"group{g}_{sn}"
            out_tsv = os.path.join(chim_dir, f"{prefix}_{tool}_blastn.tsv")
            rpt_tsv = os.path.join(chim_dir, f"{prefix}_{tool}_chimeric.tsv")
            if not args.force and os.path.exists(rpt_tsv) and os.path.getsize(rpt_tsv) > 0:
                skipped_done += 1
                continue
            blast_tasks.append((fpath, tool, sn, out_tsv, rpt_tsv,
                               db_name, args.chimeric_script, args.threads))

    if skipped_no_ref:
        print(f"⚠️  {skipped_no_ref} 个样本无对应 Group 参考")
    print(f"📋 任务统计: {len(blast_tasks)} 个新任务, {skipped_done} 个已完成")

    if not blast_tasks:
        print("⏭️  所有任务已完成")
    else:
        print(f"🔬 并行 BLASTn: {len(blast_tasks)} 任务, j={args.jobs}")
        done, ok, err = 0, 0, 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futures = {ex.submit(run_blastn_and_detect, t): t[2] for t in blast_tasks}
            pbar = tqdm(total=len(blast_tasks), desc="BLASTn+Chimeric", unit="task")
            for fut in as_completed(futures):
                msg = fut.result()
                if msg.startswith("✅"):
                    ok += 1
                else:
                    err += 1
                done += 1
                pbar.set_postfix_str(f"OK={ok} err={err}")
                pbar.update(1)
            pbar.close()
        print(f"完成: {ok} 成功, {err} 失败")

    # ---- 4. 解析汇总 ----
    print("📊 解析嵌合报告...")
    rows = []
    for s in samples:
        sn = s['base_name']
        info = {'Mutation_Rate': s.get('mutation_rate', '0.0'),
                'Depth': s.get('depth', 0), 'Rep': s.get('rep', 'rep1'),
                'Group': s.get('group', 0),
                'Background_Ratio': s.get('background_ratio', None)}
        files = get_assembly_files(s['sample_dir'], sn, args.mode)
        for tool, fpath in files:
            g = s.get('group', 0)
            prefix = f"group{g}_{sn}"
            rpt_tsv = os.path.join(chim_dir, f"{prefix}_{tool}_chimeric.tsv")
            total = count_contigs_in_fasta(fpath)
            n_chim = count_chimeric(rpt_tsv)
            rate = (n_chim / total * 100) if total > 0 else 0.0
            row = info.copy()
            row.update({
                'Tool': tool,
                'Total_Contigs': total,
                'Chimeric_Count': n_chim,
                'Chimeric_Rate_pct': round(rate, 2),
            })
            rows.append(row)

    if not rows:
        print("⚠️  未解析到任何数据")
        return 1

    df = pd.DataFrame(rows)
    out_file = os.path.join(args.out_dir, "chimeric_summary.tsv")
    df.to_csv(out_file, sep='\t', index=False)
    print(f"✅ chimeric_summary.tsv: {len(df)} 行, {df['Tool'].nunique()} 个工具")
    for tool in sorted(df['Tool'].unique()):
        td = df[df['Tool'] == tool]
        print(f"   {tool}: avg_chimeric={td['Chimeric_Rate_pct'].mean():.1f}% "
              f"avg_contigs={td['Total_Contigs'].mean():.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
