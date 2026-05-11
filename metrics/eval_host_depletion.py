#!/usr/bin/env python3
"""
eval_host_depletion.py (v3) — 宿主过滤消融实验评估
  使用 ripgrep (rg) 高速扫描 FASTQ，多进程并行，tqdm 进度条。

核心逻辑：
  从文件名提取病毒 Accession (如 U56975.1)，
  用 rg 统计该病毒 reads 数，其余全算宿主。

用法:
  python eval_host_depletion.py \
      --d0 step5_host_free/D0_baseline/ \
      --d1 step5_host_free/D1_kraken2_only/ \
      --d2 step5_host_free/D2_hisat2_only/ \
      --d3 step5_host_free/D3_k2_hisat2/ \
      --d4 step5_host_free/D4_full/ \
      --output host_depletion_report.tsv \
      --threads 8
"""

import os, sys, re, glob, gzip, argparse, subprocess
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
import pandas as pd
from tqdm import tqdm


def extract_virus_accession(filename):
    """从 LoD 文件名提取病毒 Accession"""
    basename = os.path.basename(filename)
    m = re.search(r'LoD_Mixed_([A-Z]{1,2}_?\d{5,}\.\d{1,2})', basename)
    if m: return m.group(1)
    m = re.search(r'LoD_Mixed_([A-Za-z]{1,4}\d{5,}\.\d{1,2})', basename)
    return m.group(1) if m else None


def count_with_rg(fastq_file, virus_acc):
    """用 ripgrep -z (搜索压缩文件) 统计 FASTQ 中宿主/病毒 reads 数"""
    try:
        vc = subprocess.run(
            ["rg", "-z", "-c", f"^@{virus_acc}-", fastq_file],
            capture_output=True, text=True, timeout=120
        )
        virus_count = int(vc.stdout.strip()) if vc.stdout.strip().isdigit() else 0
    except:
        virus_count = 0
    try:
        tc = subprocess.run(
            ["rg", "-z", "-c", "^@", fastq_file],
            capture_output=True, text=True, timeout=120
        )
        total_count = int(tc.stdout.strip()) if tc.stdout.strip().isdigit() else 0
    except:
        total_count = 0
    return max(0, total_count - virus_count), virus_count


def process_file(args_tuple):
    fastq_file, virus_acc, group_name = args_tuple
    host, virus = count_with_rg(fastq_file, virus_acc)
    basename = os.path.basename(fastq_file)
    lod_match = re.search(r'_(\d+\.\d+)', basename.replace('.fastq', '').replace('.gz', ''))
    lod_factor = float(lod_match.group(1)) if lod_match else None
    total = host + virus
    return {
        'Group': group_name, 'File': basename, 'Accession': virus_acc,
        'LoD_Factor': lod_factor, 'Host_Reads': host, 'Virus_Reads': virus,
        'Total_Reads': total,
        'Virus_Pct': round(virus / total * 100, 4) if total > 0 else 0,
    }


def scan_dir(dir_path, group_name, threads):
    files = sorted(glob.glob(os.path.join(dir_path, "*_R1.fastq*")))
    tasks = []
    for f in files:
        acc = extract_virus_accession(f)
        if acc is None:
            print(f"  WARNING: Cannot extract accession from {os.path.basename(f)}, skipping")
            continue
        tasks.append((f, acc, group_name))

    print(f"[{group_name}] {len(tasks)} files, {threads} threads")
    records = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        records = list(tqdm(executor.map(process_file, tasks), total=len(tasks), desc=group_name))

    total_host = sum(r['Host_Reads'] for r in records)
    total_virus = sum(r['Virus_Reads'] for r in records)
    print(f"  Host: {total_host:,}  Virus: {total_virus:,}")
    return records


def main():
    parser = argparse.ArgumentParser(description="宿主过滤消融实验评估 (v3, rg+并行)")
    parser.add_argument('--d0', required=True)
    parser.add_argument('--d1', required=True)
    parser.add_argument('--d2', required=True)
    parser.add_argument('--d3', required=True)
    parser.add_argument('--d4', required=True)
    parser.add_argument('--outdir', default='step5_host_free_analysis', help='输出目录 (default: step5_host_free_analysis)')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 检查 ripgrep
    if not subprocess.run(["which", "rg"], capture_output=True).returncode == 0:
        print("ERROR: ripgrep (rg) not found. Install: conda install -c conda-forge ripgrep")
        sys.exit(1)

    all_records = []
    all_records.extend(scan_dir(args.d0, 'D0_Baseline', args.threads))
    all_records.extend(scan_dir(args.d1, 'D1_Kraken2_only', args.threads))
    all_records.extend(scan_dir(args.d2, 'D2_HISAT2_only', args.threads))
    all_records.extend(scan_dir(args.d3, 'D3_K2_HISAT2', args.threads))
    all_records.extend(scan_dir(args.d4, 'D4_Full', args.threads))

    detail_df = pd.DataFrame(all_records)

    # 计算每个样本的病毒保留率（以D0为基线）
    d0_virus = detail_df[detail_df['Group'] == 'D0_Baseline'].set_index('Accession')['Virus_Reads'].to_dict()
    def calc_retention(row):
        if row['Group'] == 'D0_Baseline': return 100.0
        base = d0_virus.get(row['Accession'], None)
        return round(row['Virus_Reads'] / base * 100, 2) if base and base > 0 else None
    detail_df['Virus_Retention'] = detail_df.apply(calc_retention, axis=1)

    # detail 表
    detail_path = os.path.join(args.outdir, 'host_depletion_detail.tsv')
    detail_df.to_csv(detail_path, sep='\t', index=False)

    # 汇总表
    summary = detail_df.groupby('Group').agg(Host_Total=('Host_Reads', 'sum'), Virus_Total=('Virus_Reads', 'sum')).reset_index()
    bh, bv = summary.iloc[0]['Host_Total'], summary.iloc[0]['Virus_Total']
    summary['Host_Depletion_Rate(%)'] = ((1 - summary['Host_Total'] / bh) * 100).round(2)
    summary['Virus_Retention_Rate(%)'] = ((summary['Virus_Total'] / bv) * 100).round(2)
    br = bv / (bh + bv) if (bh + bv) > 0 else 1
    summary['Virus_Enrichment_Fold'] = summary.apply(
        lambda r: round((r['Virus_Total'] / (r['Host_Total'] + r['Virus_Total'])) / br, 1)
        if (r['Host_Total'] + r['Virus_Total']) > 0 else 0, axis=1)

    summary_path = os.path.join(args.outdir, 'host_depletion_summary.tsv')
    summary.to_csv(summary_path, sep='\t', index=False)

    print("\n" + summary.to_string(index=False))
    print(f"\nDetail: {detail_path}")
    print(f"Summary: {summary_path}")
    print(f"\nPlot: python eval_plant_virus/metrics/plot_host_depletion.py --detail-tsv {detail_path} --summary-tsv {summary_path} --outdir {args.outdir}")


if __name__ == "__main__":
    main()
