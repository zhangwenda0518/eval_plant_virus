#!/usr/bin/env python3
"""
eval_host_depletion.py (all-in-one 支持版 - 极速优化)
— 自动探测 D0~D4 目录、支持 all-in-one 模式每个病毒独立统计、多进程并行。
— 引入 awk 单次流式扫描，大幅提升巨型 FASTQ 文件的解析速度。
"""

import os, sys, re, glob, argparse, subprocess
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from tqdm import tqdm


def parse_filename(basename):
    """
    提取 Accession 和 LoD_Factor
    支持:
      - LoD_Mixed_{ACCESSION}_{LoD}x_PE_[clean_]R1.fastq.gz
      - Host_Depletion_Mixed_PE_[clean_]R1.fastq.gz  (all-in-one, LoD=1000)
    """
    # 匹配: LoD_Mixed_ + Accession + _ + LoD数字 + x
    m = re.search(r'LoD_Mixed_([A-Za-z]{1,4}_?\d{4,}\.\d{1,2})_(\d+\.?\d*)x?_', basename)
    if m:
        return m.group(1), float(m.group(2))
    # 兼容旧版: LoD_Mixed_xxx_0.0005_PE_R1
    m = re.search(r'LoD_Mixed_([A-Za-z]{1,4}_?\d{4,}\.\d{1,2})_(\d+\.\d+)_', basename)
    if m:
        return m.group(1), float(m.group(2))
    # Host_Depletion_Mixed (all-in-one): 返回占位符，后续从金标准加载病毒列表
    if basename.startswith('Host_Depletion_Mixed'):
        return 'ALL_IN_ONE', 1000.0
    return None, None


def load_virus_list(gold_tsv=None):
    """从 SpikeIn_GroundTruth.tsv 加载病毒 Accession 列表"""
    search_paths = [
        gold_tsv,
        'step2_simulator/eval1_host_dep/SpikeIn_GroundTruth.tsv',
        'SpikeIn_GroundTruth.tsv',
    ]
    for p in search_paths:
        if p and os.path.exists(p):
            df = pd.read_csv(p, sep='\t')
            if 'Virus_File' in df.columns:
                viruses = [f.replace('.fasta', '').replace('.fa', '') for f in df['Virus_File']]
                print(f"[Gold] Loaded {len(viruses)} viruses from {p}")
                return viruses
    print("[Gold] Warning: No SpikeIn_GroundTruth.tsv found, using aggregate mode")
    return []


def count_reads(fastq_file, virus_acc):
    """统计指定病毒的 reads 数和总 reads 数 (传统模式)"""
    virus_count, total_count = 0, 0
    basename = os.path.basename(fastq_file)

    try:
        vc = subprocess.run(
            ["rg", "-z", "-c", f"^@{virus_acc}-", fastq_file],
            capture_output=True, text=True
        )
        if vc.stdout.strip().isdigit():
            virus_count = int(vc.stdout.strip())
    except Exception as e:
        print(f"\n[Error] 扫描病毒 Reads 失败 {basename}: {e}")

    try:
        tc_cmd = f"zcat '{fastq_file}' 2>/dev/null | wc -l"
        tc = subprocess.run(tc_cmd, shell=True, capture_output=True, text=True)
        if not tc.stdout.strip().isdigit() or int(tc.stdout.strip()) == 0:
            tc_cmd = f"gzcat '{fastq_file}' 2>/dev/null | wc -l"
            tc = subprocess.run(tc_cmd, shell=True, capture_output=True, text=True)
        if tc.stdout.strip().isdigit():
            total_count = int(tc.stdout.strip()) // 4
    except Exception as e:
        print(f"\n[Error] 统计总 Reads 失败 {basename}: {e}")

    host_count = max(0, total_count - virus_count)
    return host_count, virus_count


def count_all_reads(fastq_file):
    """All-in-one 汇总模式 (无病毒列表时的 Fallback)"""
    virus_count, host_count = 0, 0
    try:
        vc = subprocess.run(["rg", "-z", "-c", r'^@[A-Z]{1,4}_?\d{4,}\.\d{1,2}-', fastq_file], capture_output=True, text=True)
        if vc.stdout.strip().isdigit(): virus_count = int(vc.stdout.strip())
    except: pass

    try:
        hc = subprocess.run(["rg", "-z", "-c", r'^@chr', fastq_file], capture_output=True, text=True)
        if hc.stdout.strip().isdigit(): host_count = int(hc.stdout.strip())
    except: pass

    try:
        tc = subprocess.run(f"zcat '{fastq_file}' 2>/dev/null | wc -l", shell=True, capture_output=True, text=True)
        if tc.stdout.strip().isdigit():
            total = int(tc.stdout.strip()) // 4
            unclassified = max(0, total - virus_count - host_count)
            host_count += unclassified
    except: pass

    return host_count, virus_count


def count_all_in_one_fast(fastq_file):
    """
    极速模式：一次性流式解压并扫描 FASTQ 文件，在内存中完成所有前缀的 Reads 汇总。
    """
    # awk 逻辑: 只看 NR%4==1 (Header行)
    cmd = f"(zcat '{fastq_file}' 2>/dev/null || gzcat '{fastq_file}' 2>/dev/null) | awk -F'[@-]' 'NR%4==1 {{counts[$2]++}} END {{for(k in counts) print counts[k], k}}'"
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    except Exception as e:
        print(f"\n[Error] 极速扫描文件失败 {fastq_file}: {e}")
        return {}, 0
        
    counts = {}
    total_reads = 0
    for line in result.stdout.strip().split('\n'):
        if not line: continue
        parts = line.strip().split()
        if len(parts) >= 2:
            count = int(parts[0])
            prefix = parts[1]
            counts[prefix] = count
            total_reads += count

    return counts, total_reads


def process_file(args_tuple):
    fastq_file, group_name, virus_list = args_tuple
    basename = os.path.basename(fastq_file)
    virus_acc, lod_factor = parse_filename(basename)

    if not virus_acc:
        return None

    if virus_acc == 'ALL_IN_ONE' and virus_list:
        # 🚀 极速扫描替代循环
        counts_dict, total_reads = count_all_in_one_fast(fastq_file)
        
        # 🌟 修正点：先计算文件中所有已知病毒的总 Reads
        total_known_virus_reads = sum(counts_dict.get(v, 0) for v in virus_list)
        # 真正的宿主 Reads = 总 Reads - 所有病毒的 Reads (保持文件级别唯一)
        real_host_reads = max(0, total_reads - total_known_virus_reads)
        
        results = []
        for vacc in virus_list:
            virus_count = counts_dict.get(vacc, 0)
            total = max(1, total_reads)
            
            results.append({
                'Group': group_name,
                'File': basename,
                'Accession': vacc,
                'LoD_Factor': lod_factor,
                'Host_Reads': real_host_reads,  # 改用真实的全局宿主Reads
                'Virus_Reads': virus_count,
                'Total_Reads': total,
                'Virus_Pct': round(virus_count / total * 100, 4),
            })
        return results
        
    elif virus_acc == 'ALL_IN_ONE':
        # fallback: 无病毒列表时汇总统计
        host, virus = count_all_reads(fastq_file)
        total = max(1, host + virus)
        return [{
            'Group': group_name, 'File': basename,
            'Accession': virus_acc, 'LoD_Factor': lod_factor,
            'Host_Reads': host, 'Virus_Reads': virus,
            'Total_Reads': total,
            'Virus_Pct': round(virus / total * 100, 4),
        }]
        
    else:
        # 标准单独病毒 LoD 模式
        host, virus = count_reads(fastq_file, virus_acc)
        total = max(1, host + virus)
        return [{
            'Group': group_name, 'File': basename,
            'Accession': virus_acc, 'LoD_Factor': lod_factor,
            'Host_Reads': host, 'Virus_Reads': virus,
            'Total_Reads': total,
            'Virus_Pct': round(virus / total * 100, 4),
        }]


def scan_dir(dir_path, group_name, jobs, virus_list=None):
    files = sorted(glob.glob(os.path.join(dir_path, "*_R1.fastq*")))
    tasks = [(f, group_name, virus_list) for f in files]

    if not tasks:
        print(f"\n[Warning] 在 {dir_path} 中没有找到任何 FASTQ 文件！")
        return []

    print(f"[{group_name}] 正在使用 {jobs} 个并发进程扫描 {len(tasks)} 个文件...")
    records = []

    with ProcessPoolExecutor(max_workers=jobs) as executor:
        results = list(tqdm(executor.map(process_file, tasks), total=len(tasks), desc=group_name))

    for res in results:
        if res is not None:
            if isinstance(res, list):
                records.extend(res)
            else:
                records.append(res)

    total_host = sum(r['Host_Reads'] for r in records)
    total_virus = sum(r['Virus_Reads'] for r in records)
    print(f"  -> 完成! 提取到 {len(records)} 条记录. 总宿主: {total_host:,} | 总病毒: {total_virus:,}\n")
    return records


def main():
    parser = argparse.ArgumentParser(description="宿主过滤消融实验评估 (极速版)")
    parser.add_argument('-i', '--input_dir', required=True,
                        help='包含 D0~D4 文件夹的父目录 (如: step5_host_free)')
    parser.add_argument('--outdir', default='step5_host_free_analysis', help='输出目录')
    parser.add_argument('-j', '--jobs', type=int, default=8, help='并行运行的进程数')
    parser.add_argument('--gold', default=None,
                        help='金标准文件路径 (默认自动搜索 SpikeIn_GroundTruth.tsv)')
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"[Error] 输入目录不存在: {args.input_dir}")
        sys.exit(1)

    dir_map = {'D0': None, 'D1': None, 'D2': None, 'D3': None, 'D4': None}
    for folder in os.listdir(args.input_dir):
        folder_path = os.path.join(args.input_dir, folder)
        if not os.path.isdir(folder_path):
            continue
        prefix = folder[:2].upper()
        if prefix in dir_map:
            dir_map[prefix] = folder_path

    missing = [k for k, v in dir_map.items() if v is None]
    if missing:
        print(f"[Error] 探测失败！在 {args.input_dir} 中找不到以下前缀的文件夹: {', '.join(missing)}")
        sys.exit(1)

    print("\n✅ 成功检测到目录结构:")
    for k in sorted(dir_map.keys()):
        print(f"  {k} -> {os.path.basename(dir_map[k])}")
    print("-" * 50 + "\n")

    os.makedirs(args.outdir, exist_ok=True)
    virus_list = load_virus_list(args.gold)

    all_records = []
    scan_specs = [
        ('D0', 'D0_Baseline'),
        ('D1', 'D1_Kraken2_only'),
        ('D2', 'D2_HISAT2_only'),
        ('D3', 'D3_K2_HISAT2'),
        ('D4', 'D4_Full'),
    ]
    for key, label in scan_specs:
        if key in dir_map:
            all_records.extend(scan_dir(dir_map[key], label, args.jobs, virus_list))

    if not all_records:
        print("[Error] 没有扫描到任何数据！")
        sys.exit(1)

    detail_df = pd.DataFrame(all_records)

    # 计算指标
    d0_df = detail_df[detail_df['Group'] == 'D0_Baseline']
    d0_by_key = {(r['Accession'], r['LoD_Factor']): r['Virus_Reads'] for _, r in d0_df.iterrows()}

    def calc_retention(row):
        if row['Group'] == 'D0_Baseline': return 100.0
        base = d0_by_key.get((row['Accession'], row['LoD_Factor']), None)
        return round(row['Virus_Reads'] / base * 100, 2) if base and base > 0 else None

    detail_df['Virus_Retention'] = detail_df.apply(calc_retention, axis=1)

    d0_pct = {(r['Accession'], r['LoD_Factor']): r['Virus_Pct'] for _, r in d0_df.iterrows()}
    detail_df['Enrichment_Fold'] = detail_df.apply(
        lambda r: round(r['Virus_Pct'] / d0_pct.get((r['Accession'], r['LoD_Factor']), 0.01), 2)
        if r['Group'] != 'D0_Baseline' else 1.0, axis=1
    )

    detail_df['Host_Removal_Pct'] = detail_df.apply(
        lambda r: round((1 - r['Host_Reads'] / r['Total_Reads']) * 100, 2)
        if r['Total_Reads'] > 0 else 0, axis=1
    )

    detail_path = os.path.join(args.outdir, 'host_depletion_detail.tsv')
    detail_df.to_csv(detail_path, sep='\t', index=False)
    print(f"🎉 详细结果已保存至: {detail_path} ({len(detail_df)} 行)")

    summary = detail_df.groupby('Group').agg(
        Samples=('File', 'nunique'),
        Viruses=('Accession', 'nunique'),
        Mean_Virus_Retention=('Virus_Retention', 'mean'),
        Median_Virus_Retention=('Virus_Retention', 'median'),
        Mean_Host_Removal_Pct=('Host_Removal_Pct', 'mean'),
        Mean_Enrichment_Fold=('Enrichment_Fold', 'mean'),
        Total_Host_Reads=('Host_Reads', 'sum'),
        Total_Virus_Reads=('Virus_Reads', 'sum'),
    ).reset_index()
    
    summary_path = os.path.join(args.outdir, 'host_depletion_summary.tsv')
    summary.to_csv(summary_path, sep='\t', index=False)
    print(f"📊 汇总结果已保存至: {summary_path}")

if __name__ == "__main__":
    main()
