#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: detect_chimeric_contigs.py
功能描述: 评估病毒组装结果中的嵌合 Contig (跨不同参考基因组的错误拼接)
依赖库: pandas

用法:
  python detect_chimeric_contigs.py -i contigs_to_ref.blastn.tsv -o chimeric_report.tsv
"""

import pandas as pd
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="根据 BLASTn 结果检测嵌合 Contig")
    parser.add_argument("-i", "--input", required=True, help="BLASTn 输出结果 (outfmt 6)")
    parser.add_argument("-m", "--min_len", type=int, default=200, help="最小比对长度阈值 (默认: 200bp)")
    parser.add_argument("-p", "--min_id", type=float, default=90.0, help="最小一致性百分比 (默认: 90.0)")
    parser.add_argument("-o", "--max_overlap", type=int, default=100, help="允许的最大重叠区域，防止保守域假阳性 (默认: 100bp)")
    parser.add_argument("--out", default="chimeric_report.tsv", help="输出的嵌合事件报告")
    return parser.parse_args()

def check_chimeric(df_group, max_overlap):
    """
    检查单个 Contig 是否为嵌合体
    df_group: 同一个 Contig 的所有比对记录
    返回: bool (是否嵌合), str (嵌合细节描述)
    """
    unique_refs = df_group['sseqid'].unique()

    if len(unique_refs) < 2:
        return False, ""

    spans = {}
    for ref in unique_refs:
        ref_hits = df_group[df_group['sseqid'] == ref]
        min_start = ref_hits[['qstart', 'qend']].min(axis=1).min()
        max_end   = ref_hits[['qstart', 'qend']].max(axis=1).max()
        spans[ref] = (min_start, max_end)

    refs_list = list(spans.keys())
    for i in range(len(refs_list)):
        for j in range(i + 1, len(refs_list)):
            refA = refs_list[i]
            refB = refs_list[j]

            startA, endA = spans[refA]
            startB, endB = spans[refB]

            overlap = max(0, min(endA, endB) - max(startA, startB))

            if overlap <= max_overlap:
                detail = f"{refA}[{startA}-{endA}] AND {refB}[{startB}-{endB}] (Overlap: {overlap}bp)"
                return True, detail

    return False, ""

def main():
    args = parse_args()

    columns = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    try:
        df = pd.read_csv(args.input, sep='\t', names=columns)
    except FileNotFoundError:
        print(f"Error: 找不到文件 {args.input}", file=sys.stderr)
        sys.exit(1)

    df_filtered = df[(df['length'] >= args.min_len) & (df['pident'] >= args.min_id)].copy()

    total_contigs = df_filtered['qseqid'].nunique()
    if total_contigs == 0:
        print("警告: 过滤后没有留下任何有效的 Contig 比对结果。请检查输入或降低过滤阈值。")
        sys.exit(0)

    chimeric_results = []

    for contig_id, group in df_filtered.groupby('qseqid'):
        is_chimera, detail = check_chimeric(group, args.max_overlap)
        if is_chimera:
            chimeric_results.append({
                'Contig_ID': contig_id,
                'Chimeric_Detail': detail
            })

    num_chimeras = len(chimeric_results)
    chimeric_rate = (num_chimeras / total_contigs) * 100

    print("="*40)
    print("        嵌合 Contig 检测报告        ")
    print("="*40)
    print(f"有效评估的 Contigs 总数 : {total_contigs}")
    print(f"检测到的嵌合 Contigs 数   : {num_chimeras}")
    print(f"嵌合率 (Chimeric Rate)  : {chimeric_rate:.2f} %")
    print("="*40)

    if num_chimeras > 0:
        out_df = pd.DataFrame(chimeric_results)
        out_df.to_csv(args.out, sep='\t', index=False)
        print(f"嵌合细节已保存至 : {args.out}")

if __name__ == "__main__":
    main()
