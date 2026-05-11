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
    """
    用 ripgrep 统计一个 FASTQ 中的宿主/病毒 reads 数。
    rg -c 返回匹配行数（即 reads 数），速度极快。
    """
    # 病毒 reads：read ID 以 virus_acc 开头
    try:
        # 统计病毒 reads：rg -c "^@病毒Accession" fastq_file
        virus_cmd = ["rg", "-c", rf"^@{virus_acc}[-\s/]", fastq_file]
        virus_result = subprocess.run(virus_cmd, capture_output=True, text=True, timeout=60)
        virus_count = int(virus_result.stdout.strip()) if virus_result.stdout.strip().isdigit() else 0
    except Exception:
        virus_count = 0

    # 总 reads = @ 开头的行数
    try:
        total_cmd = ["rg", "-c", r"^@", fastq_file]
        total_result = subprocess.run(total_cmd, capture_output=True, text=True, timeout=60)
        total_count = int(total_result.stdout.strip()) if total_result.stdout.strip().isdigit() else 0
    except Exception:
        total_count = 0

    host_count = max(0, total_count - virus_count)
    return host_count, virus_count


def process_file(args_tuple):
    """单个文件的处理任务（多进程）"""
    fastq_file, virus_acc = args_tuple
    return count_with_rg(fastq_file, virus_acc)


def scan_dir(dir_path, group_name, threads):
    """扫描目录，并行统计所有 R1 文件"""
    files = sorted(glob.glob(os.path.join(dir_path, "*_R1.fastq*")))

    # 构建任务列表
    tasks = []
    for f in files:
        acc = extract_virus_accession(f)
        if acc is None:
            print(f"  WARNING: Cannot extract accession from {os.path.basename(f)}, skipping")
            continue
        tasks.append((f, acc))

    print(f"[{group_name}] {len(tasks)} files, {threads} threads")

    total_host = 0
    total_virus = 0

    with ProcessPoolExecutor(max_workers=threads) as executor:
        results = list(tqdm(
            executor.map(process_file, tasks),
            total=len(tasks),
            desc=group_name,
            unit="file"
        ))
        for h, v in results:
            total_host += h
            total_virus += v

    print(f"  Host: {total_host:,}  Virus: {total_virus:,}")
    return {'Group': group_name, 'Host_Reads': total_host, 'Virus_Reads': total_virus}


def main():
    parser = argparse.ArgumentParser(description="宿主过滤消融实验评估 (v3, rg+并行)")
    parser.add_argument('--d0', required=True)
    parser.add_argument('--d1', required=True)
    parser.add_argument('--d2', required=True)
    parser.add_argument('--d3', required=True)
    parser.add_argument('--d4', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()

    # 检查 ripgrep
    if not subprocess.run(["which", "rg"], capture_output=True).returncode == 0:
        print("ERROR: ripgrep (rg) not found. Install: conda install -c conda-forge ripgrep")
        sys.exit(1)

    results = [
        scan_dir(args.d0, 'D0_Baseline', args.threads),
        scan_dir(args.d1, 'D1_Kraken2_only', args.threads),
        scan_dir(args.d2, 'D2_HISAT2_only', args.threads),
        scan_dir(args.d3, 'D3_K2_HISAT2', args.threads),
        scan_dir(args.d4, 'D4_Full', args.threads),
    ]

    df = pd.DataFrame(results)
    base_host = df.loc[df['Group'] == 'D0_Baseline', 'Host_Reads'].values[0]
    base_virus = df.loc[df['Group'] == 'D0_Baseline', 'Virus_Reads'].values[0]

    df['Host_Depletion_Rate(%)'] = ((1 - df['Host_Reads'] / base_host) * 100).round(2)
    df['Virus_Retention_Rate(%)'] = ((df['Virus_Reads'] / base_virus) * 100).round(2)

    base_ratio = base_virus / (base_host + base_virus) if (base_host + base_virus) > 0 else 1
    def enrichment(row):
        total = row['Host_Reads'] + row['Virus_Reads']
        return round((row['Virus_Reads'] / total) / base_ratio, 1) if total > 0 else 0.0
    df['Virus_Enrichment_Fold'] = df.apply(enrichment, axis=1)

    df.to_csv(args.output, sep='\t', index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()
