#!/usr/bin/env python3
"""
计算 PR 曲线下面积 (AUPRC) 及绘制精确率-召回率曲线

用于评估三（候选病毒鉴定策略比较）的不平衡数据集评估。
AUPRC 比单一 F1 值更具说服力，适用于正负样本不平衡的场景。

用法:
  python calc_auprc.py --predictions results.tsv --labels sequence_labels.tsv --out auprc_results/
"""

import argparse, os
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score, auc


def load_labels(labels_tsv):
    """加载真值标签，返回 {seq_id: 1(positive)/0(negative)}"""
    df = pd.read_csv(labels_tsv, sep='\t')
    labels = {}
    for _, row in df.iterrows():
        labels[row['seq_id']] = 1 if row.get('type', '') == 'positive' else 0
    return labels


def load_predictions(pred_tsv, score_col='bitscore'):
    """
    加载鉴定结果。
    期望格式: seq_id \t [各种列] \t bitscore (或 e-value)
    如果有 score 则用 score，否则用 e-value 的 -log10 转换
    """
    df = pd.read_csv(pred_tsv, sep='\t')
    preds = {}
    for _, row in df.iterrows():
        seq_id = row.get('seq_id') or row.get('qseqid') or row.get('contig_id')
        if not seq_id:
            continue
        if score_col in df.columns:
            score = float(row[score_col])
        elif 'evalue' in df.columns or 'e-value' in df.columns:
            ev = float(row.get('evalue', row.get('e-value', 1)))
            score = -np.log10(max(ev, 1e-300))
        else:
            score = 1.0  # 无分数时默认为1
        if seq_id not in preds or score > preds[seq_id]:
            preds[seq_id] = score
    return preds


def compute_auprc(labels, preds):
    """计算 AUPRC"""
    y_true = []
    y_score = []

    for seq_id, label in labels.items():
        y_true.append(label)
        y_score.append(preds.get(seq_id, 0.0))

    ap = average_precision_score(y_true, y_score)
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    auc_pr = auc(recall, precision)
    return ap, auc_pr, precision, recall


def plot_pr_curve(precision, recall, ap, output_path, label="Method"):
    """绘制 PR 曲线"""
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.plot(recall, precision, linewidth=2, label=f'{label} (AUPRC={ap:.3f})')
    ax.set_xlabel('Recall', fontsize=14)
    ax.set_ylabel('Precision', fontsize=14)
    ax.set_title('Precision-Recall Curve', fontsize=16)
    ax.legend(loc='lower left', fontsize=12)
    ax.set_xlim([0.0, 1.05])
    ax.set_ylim([0.0, 1.05])
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  PR curve saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="计算 AUPRC 并绘制 PR 曲线")
    parser.add_argument("--predictions", required=True, help="鉴定结果 TSV (每行一个命中)")
    parser.add_argument("--labels", required=True, help="金标准标签 TSV")
    parser.add_argument("--score-col", default="bitscore", help="评分列名 (default: bitscore)")
    parser.add_argument("--method-label", default="Method", help="图例标签")
    parser.add_argument("--outdir", required=True, help="输出目录")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    labels = load_labels(args.labels)
    preds = load_predictions(args.predictions, args.score_col)

    ap, auc_pr, precision, recall = compute_auprc(labels, preds)

    print(f"AUPRC (average_precision): {ap:.4f}")
    print(f"AUC-PR: {auc_pr:.4f}")

    plot_pr_curve(precision, recall, ap, os.path.join(args.outdir, "pr_curve.png"), args.method_label)

    # 写结果
    with open(os.path.join(args.outdir, "auprc_summary.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"AUPRC\t{ap:.6f}\n")
        f.write(f"AUC_PR\t{auc_pr:.6f}\n")
        f.write(f"n_positive\t{sum(labels.values())}\n")
        f.write(f"n_total\t{len(labels)}\n")
    print(f"  Summary saved: {os.path.join(args.outdir, 'auprc_summary.tsv')}")


if __name__ == "__main__":
    main()
