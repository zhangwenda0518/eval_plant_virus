#!/usr/bin/env python3
"""评估三：候选病毒鉴定策略比较 — 绘图"""

import argparse, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from plot_utils import set_style, get_color, save_default

def plot_pr_curves(pr_data, outpath):
    """图7: PR曲线 + AUPRC"""
    set_style()
    fig, ax = plt.subplots(figsize=(8, 7))

    for label, (precision, recall, auprc) in pr_data.items():
        ax.plot(recall, precision, linewidth=2, label=f'{label} (AUPRC={auprc:.3f})',
                color=get_color(label))

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Curves for Identification Strategies')
    ax.legend(loc='lower left')
    ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    save_default(fig, outpath)

def plot_contribution_upset(contrib_data, outpath):
    """图8: 搜索原理贡献度柱状图（简化版UpSet替代）"""
    set_style()
    fig, ax = plt.subplots(figsize=(9, 6))

    configs = list(contrib_data.keys())
    counts = list(contrib_data.values())
    colors = ['#4C72B0'] * len(configs)
    if 'P5' in configs: colors[configs.index('P5')] = '#C44E52'

    bars = ax.bar(configs, counts, color=colors, edgecolor='white')
    ax.set_ylabel('Unique Positive Sequences Detected')
    ax.set_xlabel('Search Configuration')
    ax.set_title('Contribution of Each Search Principle')
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                str(c), ha='center')
    save_default(fig, outpath)

def plot_adversarial_threshold_heatmap(threshold_data, outpath):
    """图9: 对抗策略阈值扫描热图"""
    import seaborn as sns
    set_style()
    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(threshold_data, annot=True, fmt='.3f', cmap='YlOrRd',
                ax=ax, cbar_kws={'label': 'F1 Score'})
    ax.set_title('Adversarial Strategy: Threshold × Filter Sensitivity')
    ax.set_xlabel('Filter Parameter')
    ax.set_ylabel('B_viral / B_NR Ratio')
    save_default(fig, outpath)

def main():
    parser = argparse.ArgumentParser('评估三绘图')
    parser.add_argument('--pr-tsv', help='PR曲线数据 (method,precision,recall,auprc)')
    parser.add_argument('--contrib-tsv', help='贡献度数据 (config,count)')
    parser.add_argument('--threshold-tsv', help='阈值扫描矩阵 TSV')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.pr_tsv and os.path.exists(args.pr_tsv):
        prdf = pd.read_csv(args.pr_tsv, sep='\t')
        pr_data = {}
        for method in prdf['method'].unique():
            md = prdf[prdf['method']==method]
            pr_data[method] = (md['precision'].values, md['recall'].values, md['auprc'].values[0])
        plot_pr_curves(pr_data, os.path.join(args.outdir, 'Fig7_PR_Curves'))

    if args.contrib_tsv and os.path.exists(args.contrib_tsv):
        cdf = pd.read_csv(args.contrib_tsv, sep='\t')
        contrib = dict(zip(cdf.iloc[:,0], cdf.iloc[:,1].astype(int)))
        plot_contribution_upset(contrib, os.path.join(args.outdir, 'Fig8_Contribution'))

    if args.threshold_tsv and os.path.exists(args.threshold_tsv):
        tdf = pd.read_csv(args.threshold_tsv, sep='\t', index_col=0)
        plot_adversarial_threshold_heatmap(tdf, os.path.join(args.outdir, 'Fig9_Threshold_Heatmap'))

if __name__ == '__main__':
    main()
