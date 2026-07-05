#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Virus Assembly Benchmark — Phase 1: MetaQUAST 质量评估
================================================================================

功能:
  1. 自动扫描所有 sample 目录，识别组装文件
  2. 并行运行 MetaQUAST (7 种或 4 种组装方法)
  3. 解析 transposed_report.tsv → metaquast_summary.tsv

断点续跑: 检查 benchmark_results/metaquast/{sample}/metaquast.log 是否存在

用法:
  python scripts/benchmark_mq.py -i step6_assemblies -r ref.fasta -o benchmark_results --mode 7 -j 4 -t 4
  python scripts/benchmark_mq.py -i step6_assemblies -r ref.fasta -o benchmark_results --mode 7 --skip-run

输出:
  - benchmark_results/metaquast/{sample}/          ← MetaQUAST 原始输出
  - benchmark_results/metaquast_summary.tsv
"""

import os, sys, argparse, subprocess, shutil, re
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm

# 导入公共模块
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import (
    discover_samples, get_assembly_files, parse_dir_name, METAQUAST_METRICS
)


def run_single_metaquast(args_for_mq):
    """运行单个 sample 的 MetaQUAST，供并行调度"""
    sample_dir, base_name, out_dir, ref, min_contig, threads, mode, force, ref_dir = args_for_mq
    mq_log = os.path.join(out_dir, "metaquast.log")
    # 检查是否真的完成: metaquast.log + per-reference 报告存在
    done_marker = os.path.join(out_dir, "runs_per_reference")
    if not force and os.path.exists(mq_log) and os.path.isdir(done_marker) and os.listdir(done_marker):
        return f"⏭️ {base_name}"

    files = get_assembly_files(sample_dir, base_name, mode)
    if len(files) == 0:
        return f"❌ {base_name}: no assemblies"

    labels = ",".join(t for t, _ in files)
    paths = [p for _, p in files]
    os.makedirs(out_dir, exist_ok=True)
    # 解析 metaquast 完整路径 (在子进程中 PATH 可能不同，优先显式查找)
    metaquast_exe = (shutil.which("metaquast.py") or shutil.which("metaquast")
                     or os.path.expanduser("~/mambaforge/bin/metaquast.py")
                     or "metaquast.py")
    # 用参考目录触发 per-reference 分箱; 若无 ref_dir 则用合并 FASTA
    cmd = [
        metaquast_exe, "-o", out_dir, "-r", ref_dir if ref_dir else ref,
        "--unique-mapping",
        "--min-contig", str(min_contig), "-t", str(threads),
        "-l", labels
    ] + paths
    try:
        with open(mq_log, 'w') as lf:
            subprocess.run(cmd, check=True, stdout=lf, stderr=subprocess.STDOUT, timeout=86400)
        return f"✅ {base_name}"
    except subprocess.TimeoutExpired:
        return f"❌ {base_name}: timeout"
    except subprocess.CalledProcessError as e:
        return f"❌ {base_name}: exit={e.returncode}"
    except FileNotFoundError as e:
        return f"❌ {base_name}: metaquast not found in PATH"
    except Exception as e:
        return f"❌ {base_name}: {e}"


def parse_metaquast_reports(out_base):
    """解析所有 MetaQUAST per-reference reports → DataFrame"""
    import glob as gmod
    mq_base = os.path.join(out_base, "metaquast")

    # 两种路径模式:
    #   A: metaquast/group{g}/{sample}/runs_per_reference/{virus}/transposed_report.tsv
    #   B: metaquast/{sample}/runs_per_reference/{virus}/transposed_report.tsv
    pattern_a = os.path.join(mq_base, "group*", "*", "runs_per_reference", "*", "transposed_report.tsv")
    pattern_b = os.path.join(mq_base, "*", "runs_per_reference", "*", "transposed_report.tsv")

    rows = []
    for path in gmod.glob(pattern_a) + gmod.glob(pattern_b):
        try:
            parts = os.path.normpath(path).split(os.sep)
            virus_dir = parts[-2]
            if virus_dir == "not_aligned":
                continue
            virus_name = virus_dir.split(".")[0] if "." in virus_dir else virus_dir

            # 判断是否有 group 层级
            if parts[-5].startswith("group"):
                base_name  = parts[-4]
                group_dir  = parts[-5]
                g_match = re.search(r'group(\d+)', group_dir)
                group_id = int(g_match.group(1)) if g_match else 0
            else:
                base_name  = parts[-4]
                group_id = 0

            info = parse_dir_name(base_name)
            info.setdefault('mutation_rate', '0.0')
            info.setdefault('depth', 0)
            info.setdefault('rep', 'rep1')
            info.setdefault('background_ratio', None)

            df = pd.read_csv(path, sep='\t')
            df['Tool'] = df['Assembly']
            df['Virus'] = virus_name
            df['Mutation_Rate'] = info['mutation_rate']
            df['Depth'] = info['depth']
            df['Rep'] = info['rep']
            df['Group'] = group_id
            df['Background_Ratio'] = info['background_ratio']
            rows.append(df)
        except Exception:
            continue
    if not rows:
        return pd.DataFrame()
    df_all = pd.concat(rows, ignore_index=True)
    keep_cols = ['Tool', 'Mutation_Rate', 'Depth', 'Rep', 'Group', 'Virus', 'Background_Ratio']
    keep_cols += [c for c in METAQUAST_METRICS if c in df_all.columns]
    return df_all[[c for c in keep_cols if c in df_all.columns]]


def _split_fasta_by_group(fasta_file, out_base):
    """
    将合并的参考 FASTA 按 virus accession 中的 Group 信息拆分为 group_1/.../group_5 子目录。
    50 个 virus 的 accession 已知归属: 每 10 个一组。
    如果无法从 accession 推断 Group，则用中位数均分。
    """
    import os
    os.makedirs(out_base, exist_ok=True)

    # 读取所有序列
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # 处理 QIIME/UCSC 格式的 ID: >accession description
                full_id = line.lstrip('>').strip()
                current_id = full_id.split()[0] if full_id else "unknown"
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    if not sequences:
        return {}

    # 尝试从 step1_eval_viruses/group_N/ 推断 accession → group 映射
    # 注意: 工作目录不确定，需要尝试多个可能路径
    accession_group = {}
    possible_parents = [
        os.path.dirname(out_base),                          # benchmark_results/
        os.path.join(os.path.dirname(out_base), ".."),      # benchmark_results/../
        os.getcwd(),                                        # 当前目录
        os.path.expanduser("~/snap/simulator-data"),        # 显式指定
    ]
    for parent in possible_parents:
        for g in range(1, 6):
            gdir = os.path.join(parent, "step1_eval_viruses", f"group_{g}")
            if os.path.isdir(gdir):
                for fname in os.listdir(gdir):
                    if fname.endswith('.fasta'):
                        acc = fname.replace('.fasta', '')
                        accession_group[acc] = g
        if accession_group:
            break

    # 未匹配的按顺序均分到 5 组
    unassigned = [acc for acc in sorted(sequences) if acc not in accession_group]
    if unassigned:
        per_group = max(1, len(unassigned) // 5)
        for i, acc in enumerate(unassigned):
            g = min(5, (i // per_group) + 1)
            accession_group[acc] = g

    # 写入各 group 子目录
    group_seqs = {g: {} for g in range(1, 6)}
    for acc, seq in sequences.items():
        g = accession_group.get(acc, 1)
        group_seqs.setdefault(g, {})[acc] = seq

    group_ref_map = {}
    for g, seqs in group_seqs.items():
        if not seqs:
            continue
        gdir = os.path.join(out_base, f"group_{g}")
        os.makedirs(gdir, exist_ok=True)
        for acc, seq in seqs.items():
            with open(os.path.join(gdir, f"{acc}.fasta"), 'w') as f:
                f.write(f">{acc}\n{seq}\n")
        group_ref_map[g] = gdir

    return group_ref_map

def main():
    parser = argparse.ArgumentParser(description="🧬 Phase 1: MetaQUAST 质量评估")
    parser.add_argument("-i", "--input-dir", default="step6_assemblies")
    parser.add_argument("-r", "--ref", required=True,
                        help="参考病毒: 可以是 (1) 包含 group_1/... 子目录的父目录, "
                             "或 (2) 合并的 FASTA 文件，自动按 Group 拆分")
    parser.add_argument("-o", "--out-dir", default="benchmark_results")
    parser.add_argument("--mode", choices=['4', '6', '7', '10'], default='6', help="4/6/7/10 种组装方法")
    parser.add_argument("-j", "--jobs", type=int, default=4, help="并行样本数")
    parser.add_argument("-t", "--threads", type=int, default=4, help="每个 MetaQUAST 线程数")
    parser.add_argument("-m", "--min-contig", type=int, default=200, help="最小 contig 长度")
    parser.add_argument("--skip-run", action="store_true", help="跳过 MetaQUAST，只解析已有结果")
    parser.add_argument("--force", action="store_true", help="强制重新运行")

    # 🌟 新增的参数
    parser.add_argument("--split-ref", choices=['auto', 'always', 'never'], default='auto',
                        help="处理合并的 FASTA 时是否按 Group 拆分 (默认 auto: 根据样本自动判断)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ---- 1. 提前扫描样本 (用于自动判断是否需要拆分) ----
    samples = discover_samples(args.input_dir)
    if not samples:
        print("❌ 未发现任何样本目录！")
        return 1
    print(f"🔍 发现 {len(samples)} 个样本")

    # 智能判断：如果解析出的样本中，有任何一个的 group 不是 0，说明存在 group 结构
    has_groups = any(s.get('group', 0) != 0 for s in samples)

    # ---- 2. 检查依赖 & 准备参考目录 ----
    if not args.skip_run:
        metaquast_exe = (shutil.which("metaquast.py") or shutil.which("metaquast")
                         or os.path.expanduser("~/mambaforge/bin/metaquast.py")
                         or "metaquast.py")
        print(f"✅ metaquast: {metaquast_exe}")

    ref_path = os.path.abspath(args.ref)
    group_ref_map = {}

    if os.path.isfile(ref_path):
        # 🌟 决定拆分策略
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

            # 高效流式读取：把合并的 fasta 拆成一个个独立的病毒文件，全放进 flat_ref_dir
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

            # 将默认的 group 0 指向这个包含了所有独立病毒的目录
            group_ref_map = {0: flat_ref_dir}
            print(f"📁 独立病毒拆分完成，已全部放入全局参考目录。")

    elif os.path.isdir(ref_path):
        print(f"📁 输入为参考目录: {ref_path}")
        for gdir in sorted(os.listdir(ref_path)):
            gpath = os.path.join(ref_path, gdir)
            if os.path.isdir(gpath) and gdir.startswith("group_"):
                try:
                    gnum = int(gdir.split("_")[1])
                    group_ref_map[gnum] = gpath
                except ValueError:
                    continue
        if not group_ref_map:
            # 无 group_* 子目录 → 找 FASTA 文件并拆分为独立病毒文件
            fasta_files = [f for f in os.listdir(ref_path) if f.endswith(('.fasta', '.fa', '.fna'))]
            if fasta_files:
                print(f"📦 未找到 group_* 子目录，检测到 {len(fasta_files)} 个 FASTA，正在按病毒拆分...")
                flat_ref_dir = os.path.join(args.out_dir, "_ref_split_flat")
                os.makedirs(flat_ref_dir, exist_ok=True)
                for ff in fasta_files:
                    ff_path = os.path.join(ref_path, ff)
                    with open(ff_path, 'r') as fh:
                        cid, cseq = None, []
                        for line in fh:
                            if line.startswith('>'):
                                if cid:
                                    with open(os.path.join(flat_ref_dir, f"{cid}.fasta"), 'w') as out:
                                        out.write(f">{cid}\n{''.join(cseq)}")
                                cid = line.lstrip('>').strip().split()[0]
                                cseq = []
                            else:
                                cseq.append(line)
                        if cid:
                            with open(os.path.join(flat_ref_dir, f"{cid}.fasta"), 'w') as out:
                                out.write(f">{cid}\n{''.join(cseq)}")
                group_ref_map = {0: flat_ref_dir}
                print(f"📁 拆分完成: {len(os.listdir(flat_ref_dir))} 个独立病毒, 统一放入 Group 0")
            else:
                group_ref_map = {0: ref_path}
        print(f"📁 参考目录: {len(group_ref_map)} 个 group")
    else:
        print(f"❌ 参考路径不存在: {ref_path}")
        return 1

    # ---- 3. 并行跑 MetaQUAST ----
    if not args.skip_run:
        mq_base = os.path.join(args.out_dir, "metaquast")
        tasks = []
        for s in samples:
            sn = s['base_name']
            g = s.get('group', 0)

            # 🌟 统一获取参考目录 (刚才我们已经把 group 0 指向了扁平拆分目录)
            ref_dir = group_ref_map.get(g)
            if ref_dir is None:
                print(f"⚠️  {sn}: 找不到 group {g} 的参考目录，跳过")
                continue

            if g == 0:
                group_out = os.path.join(mq_base, sn)
            else:
                group_out = os.path.join(mq_base, f"group{g}", sn)

            tasks.append((s['sample_dir'], sn, group_out,
                          ref_path, args.min_contig, args.threads, args.mode, args.force, ref_dir))
        print(f"🚀 并行跑 MetaQUAST: {len(tasks)} 个样本, j={args.jobs}, t={args.threads}")
        done, ok, skipped = 0, 0, 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futures = {ex.submit(run_single_metaquast, t): t[1] for t in tasks}
            pbar = tqdm(total=len(tasks), desc="MetaQUAST", unit="sample")
            for fut in as_completed(futures):
                msg = fut.result()
                if msg.startswith("✅"):
                    ok += 1
                elif msg.startswith("⏭️"):
                    skipped += 1
                done += 1
                pbar.set_postfix_str(f"OK={ok} skip={skipped} fail={done-ok-skipped}")
                pbar.update(1)
            pbar.close()
        print(f"MetaQUAST: {ok} 成功, {skipped} 跳过, {done-ok-skipped} 失败")

    # ---- 3. 解析报告 ----
    print("📊 解析 MetaQUAST 报告...")
    df = parse_metaquast_reports(args.out_dir)
    if df.empty:
        print("⚠️  未解析到任何 MetaQUAST 数据，请确认 MetaQUAST 是否已完成。")
        return 1

    out_file = os.path.join(args.out_dir, "metaquast_summary.tsv")
    df.to_csv(out_file, sep='\t', index=False)
    print(f"✅ metaquast_summary.tsv: {len(df)} 行, {df['Virus'].nunique()} 个病毒, {df['Tool'].nunique()} 个工具")
    return 0


if __name__ == "__main__":
    sys.exit(main())
