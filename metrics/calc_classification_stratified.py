#!/usr/bin/env python3
"""
病毒分类评估分层报告 — 按"已知种" vs "仅有属/科级参考" 分开计算准确率

用于评估四（病毒分类方法比较）的新发病毒模拟分析。

用法:
  python calc_classification_stratified.py \
      --predictions integrated_results.tsv \
      --meta test_metadata.tsv \
      --out stratified_report/
"""

import argparse, os
import pandas as pd
from collections import defaultdict


def load_predictions(pred_tsv):
    """加载整合后的分类预测结果"""
    df = pd.read_csv(pred_tsv, sep='\t')
    return df


def load_metadata(meta_tsv):
    """加载测试元数据（含真值覆盖度等）"""
    return pd.read_csv(meta_tsv, sep='\t')


def classify_stratum(meta_row):
    """
    根据覆盖度分层：
    - 100% → Known (近乎全长，分类信息完整)
    - 80% → Near-full
    - 60% → Partial
    - 40% → Fragment
    - 20% → Highly-fragmented
    """
    cov = meta_row.get('coverage_pct', 100)
    if cov >= 100:
        return "full_length"
    elif cov >= 80:
        return "near_full"
    elif cov >= 60:
        return "partial"
    elif cov >= 40:
        return "fragment"
    else:
        return "highly_fragmented"


def compute_stratified_accuracy(pred_df, meta_df, level='family'):
    """按覆盖度分层计算准确率"""
    # 建立 seq_id → 预测分类 的映射
    pred_map = {}
    for _, row in pred_df.iterrows():
        sid = row.get('contig_id') or row.get('seq_id') or row.get('id')
        if sid and level in row:
            pred_map[sid] = str(row[level]).strip()

    results = defaultdict(lambda: {"correct": 0, "total": 0, "unassigned": 0})

    for _, row in meta_df.iterrows():
        sid = row.get('seq_id') or row.get('source_accession') or row.get('id')
        stratum = classify_stratum(row)
        true_val = str(row.get(level, '')).strip()

        pred_val = pred_map.get(sid, None)

        results[stratum]["total"] += 1
        if pred_val is None or pred_val in ('', 'Unassigned', 'none', 'NA'):
            results[stratum]["unassigned"] += 1
        elif pred_val == true_val:
            results[stratum]["correct"] += 1

    # 打印
    print(f"\n{'='*50}")
    print(f"  分类层级: {level.upper()} (按覆盖度分层)")
    print(f"{'='*50}")
    rows = []
    for stratum in ['full_length', 'near_full', 'partial', 'fragment', 'highly_fragmented']:
        r = results[stratum]
        if r['total'] == 0:
            continue
        acc = r['correct'] / r['total'] * 100 if r['total'] > 0 else 0
        print(f"  {stratum:20s} | Acc: {acc:5.1f}% ({r['correct']}/{r['total']}) "
              f"| Unassigned: {r['unassigned']} ({r['unassigned']/r['total']*100:.1f}%)")
        rows.append({
            'level': level, 'stratum': stratum,
            'accuracy': round(acc, 2), 'correct': r['correct'],
            'total': r['total'], 'unassigned': r['unassigned'],
        })

    # 总体
    total_correct = sum(r['correct'] for r in results.values())
    total_all = sum(r['total'] for r in results.values())
    total_unassigned = sum(r['unassigned'] for r in results.values())
    overall = total_correct / total_all * 100 if total_all > 0 else 0
    print(f"  {'OVERALL':20s} | Acc: {overall:5.1f}% ({total_correct}/{total_all}) "
          f"| Unassigned: {total_unassigned}")

    return rows


def main():
    parser = argparse.ArgumentParser(description="病毒分类评估分层报告")
    parser.add_argument("--predictions", required=True, help="整合后分类预测 TSV")
    parser.add_argument("--meta", required=True, help="测试元数据 TSV (含 coverage_pct)")
    parser.add_argument("--outdir", required=True, help="输出目录")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    pred_df = load_predictions(args.predictions)
    meta_df = load_metadata(args.meta)

    all_rows = []
    for level in ['family', 'genus', 'species']:
        rows = compute_stratified_accuracy(pred_df, meta_df, level)
        all_rows.extend(rows)

    out_df = pd.DataFrame(all_rows)
    out_path = os.path.join(args.outdir, "stratified_accuracy.tsv")
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"\nFull report saved: {out_path}")


if __name__ == "__main__":
    main()
