#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: eval_host_depletion.py (v2)
功能: 评估宿主过滤(Host Depletion)各步骤的 Reads 保留率与去除率
核心逻辑: 从文件名提取病毒 Accession，以此精确区分宿主/病毒 reads
         （不依赖 chr/scaffold 前缀，适配任意宿主命名）

用法:
  python eval_host_depletion.py \
      --d0 step5_host_free/D0_baseline/ \
      --d1 step5_host_free/D1_kraken2_only/ \
      --d2 step5_host_free/D2_hisat2_only/ \
      --d3 step5_host_free/D3_k2_hisat2/ \
      --d4 step5_host_free/D4_full/ \
      --output eval_results/host_depletion_report.tsv
"""

import os, glob, gzip, argparse, re
import pandas as pd

def extract_virus_accession(filename):
    """
    从 LoD 文件名中提取病毒 Accession
    文件名格式: LoD_Mixed_U56975.1_0.05_PE_clean_R2.fastq.gz
    返回: U56975.1 或 None
    """
    basename = os.path.basename(filename)
    # NCBI Accession: 1-4字母 + (可选下划线) + 5+位数字 + .版本号(1-2位)
    match = re.search(r'LoD_Mixed_([A-Z]{1,2}_?\d{5,}\.\d{1,2})', basename)
    if match:
        return match.group(1)
    # 更宽泛的匹配
    match = re.search(r'LoD_Mixed_([A-Za-z]{1,4}\d{5,}\.\d{1,2})', basename)
    if match:
        return match.group(1)
    return None

def is_virus_read(read_id, virus_accession):
    """检查 read_id 是否属于指定的病毒 Accession"""
    clean_id = read_id.lstrip('@').split('-')[0].split('/')[0].split()[0]
    return clean_id == virus_accession or clean_id.startswith(virus_accession)

def count_reads(fastq_file, virus_accession):
    """统计 Fastq 文件中的宿主/病毒 reads 数"""
    host_cnt = 0
    virus_cnt = 0
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'

    try:
        with open_func(fastq_file, mode) as f:
            for line in f:
                if line.startswith('@'):
                    if is_virus_read(line, virus_accession):
                        virus_cnt += 1
                    else:
                        host_cnt += 1
                    try:
                        next(f); next(f); next(f)
                    except StopIteration:
                        break
    except Exception as e:
        print(f"Warning: {fastq_file} - {e}")
    return host_cnt, virus_cnt

def scan_dir(dir_path, group_name):
    """扫描目录下所有 R1 fastq 文件"""
    files = sorted(glob.glob(os.path.join(dir_path, "*_R1.fastq*")))
    total_host = 0
    total_virus = 0
    print(f"[{group_name}] {dir_path} -> {len(files)} files")

    for f in files:
        virus_acc = extract_virus_accession(f)
        if virus_acc is None:
            print(f"  WARNING: Cannot extract virus accession from {os.path.basename(f)}, skipping")
            continue
        h, v = count_reads(f, virus_acc)
        total_host += h
        total_virus += v
    print(f"  Host: {total_host:,}  Virus: {total_virus:,}")
    return {'Group': group_name, 'Host_Reads': total_host, 'Virus_Reads': total_virus}

def main():
    parser = argparse.ArgumentParser(description="宿主过滤消融实验效果评估 (v2)")
    parser.add_argument('--d0', required=True, help="D0_Baseline 目录")
    parser.add_argument('--d1', required=True, help="D1_Kraken2_only 目录")
    parser.add_argument('--d2', required=True, help="D2_HISAT2_only 目录")
    parser.add_argument('--d3', required=True, help="D3_K2_HISAT2 目录")
    parser.add_argument('--d4', required=True, help="D4_Full 目录")
    parser.add_argument('--output', required=True, help="输出 TSV 路径")
    args = parser.parse_args()

    results = [
        scan_dir(args.d0, 'D0_Baseline'),
        scan_dir(args.d1, 'D1_Kraken2_only'),
        scan_dir(args.d2, 'D2_HISAT2_only'),
        scan_dir(args.d3, 'D3_K2_HISAT2'),
        scan_dir(args.d4, 'D4_Full')
    ]

    df = pd.DataFrame(results)
    base_host = df.loc[df['Group'] == 'D0_Baseline', 'Host_Reads'].values[0]
    base_virus = df.loc[df['Group'] == 'D0_Baseline', 'Virus_Reads'].values[0]

    df['Host_Depletion_Rate(%)'] = ((1 - df['Host_Reads'] / base_host) * 100).round(4)
    df['Virus_Retention_Rate(%)'] = ((df['Virus_Reads'] / base_virus) * 100).round(4)

    base_ratio = base_virus / (base_host + base_virus) if (base_host + base_virus) > 0 else 1
    def enrichment(row):
        total = row['Host_Reads'] + row['Virus_Reads']
        return round((row['Virus_Reads'] / total) / base_ratio, 2) if total > 0 else 0.0
    df['Virus_Enrichment_Fold'] = df.apply(enrichment, axis=1)

    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)
    df.to_csv(args.output, sep='\t', index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nSaved: {args.output}")

if __name__ == "__main__":
    main()
