#!/usr/bin/env python3
"""
计算丰度定量准确性指标：Bray-Curtis dissimilarity, RMSE, Spearman ρ

用于评估一（已知病毒检测方法比较）的丰度估计准确性评估。
输入为预测丰度和真实丰度的配对数据。

用法:
  python calc_abundance_accuracy.py --predictions salmon_quant.tsv --gold gold_standard.tsv --out braycurtis_results/
"""

import argparse, os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import braycurtis


def load_vectors(pred_tsv, gold_csv, pred_col='TPM', gold_col='abundance'):
    """加载预测丰度和真实丰度向量（对齐到相同的病毒列表）"""
    pred_df = pd.read_csv(pred_tsv, sep='\t')
    gold_df = pd.read_csv(gold_csv, sep='\t') if gold_csv.endswith('.csv') else pd.read_csv(gold_csv, sep='\t')

    # 找到共同的病毒/分类单元
    pred_key = 'Name' if 'Name' in pred_df.columns else pred_df.columns[0]
    gold_key = 'accession' if 'accession' in gold_df.columns else gold_df.columns[0]

    common = set(pred_df[pred_key]) & set(gold_df[gold_key])
    if len(common) == 0:
        print("WARNING: No common viruses found between prediction and gold standard")

    pred_vec = []
    gold_vec = []
    for key in common:
        p = pred_df[pred_df[pred_key] == key][pred_col].values
        g = gold_df[gold_df[gold_key] == key][gold_col].values
        if len(p) > 0 and len(g) > 0:
            pred_vec.append(float(p[0]))
            gold_vec.append(float(g[0]))

    return np.array(pred_vec), np.array(gold_vec)


def compute_metrics(pred_vec, gold_vec):
    """计算三大丰度准确性指标"""
    if len(pred_vec) < 3:
        return {"error": "Too few data points (< 3)"}

    # Bray-Curtis dissimilarity (0=identical, 1=completely different)
    # 对于向量可能包含零值，需要处理
    epsilon = 1e-10
    pred_norm = pred_vec / (pred_vec.sum() + epsilon)
    gold_norm = gold_vec / (gold_vec.sum() + epsilon)
    bc = braycurtis(pred_norm, gold_norm) if pred_norm.sum() > 0 and gold_norm.sum() > 0 else 1.0

    # RMSE (均方根误差，在 log 空间以减少极端值的影响)
    rmse_log = np.sqrt(np.mean((np.log10(pred_vec + epsilon) - np.log10(gold_vec + epsilon)) ** 2))

    # L1-norm / L1误差（在归一化空间）
    l1 = np.sum(np.abs(pred_norm - gold_norm))

    # Spearman ρ
    rho, pval = stats.spearmanr(pred_vec, gold_vec)

    return {
        "bray_curtis": round(bc, 6),
        "rmse_log10": round(rmse_log, 6),
        "l1_norm": round(l1, 6),
        "spearman_rho": round(rho, 4),
        "spearman_p": round(pval, 6),
        "n_pairs": len(pred_vec),
    }


def main():
    parser = argparse.ArgumentParser(description="计算丰度定量准确性 (Bray-Curtis/RMSE/Spearman ρ)")
    parser.add_argument("--predictions", required=True, help="预测丰度 TSV")
    parser.add_argument("--gold", required=True, help="金标准丰度 TSV")
    parser.add_argument("--pred-col", default="TPM", help="预测丰度列名")
    parser.add_argument("--gold-col", default="abundance", help="金标准丰度列名")
    parser.add_argument("--outdir", required=True, help="输出目录")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    pred_vec, gold_vec = load_vectors(args.predictions, args.gold, args.pred_col, args.gold_col)
    print(f"Loaded {len(pred_vec)} paired abundance values")

    metrics = compute_metrics(pred_vec, gold_vec)

    print("=" * 40)
    print("    丰度定量准确性报告")
    print("=" * 40)
    for k, v in metrics.items():
        print(f"  {k}: {v}")

    with open(os.path.join(args.outdir, "abundance_accuracy.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        for k, v in metrics.items():
            f.write(f"{k}\t{v}\n")
    print(f"  Saved: {os.path.join(args.outdir, 'abundance_accuracy.tsv')}")


if __name__ == "__main__":
    main()
