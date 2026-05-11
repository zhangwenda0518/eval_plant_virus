#!/usr/bin/env python3
"""评估四：病毒分类方法比较 — 绘图"""

import argparse, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from plot_utils import set_style, get_color, save_default

def plot_accuracy_bar(acc_df, outpath):
    """图10: 科/属/种三级准确率对比柱状图"""
    set_style()
    methods = ['MMseqs2', 'VITAP', 'ACVirus', 'Integrated']
    levels = ['family', 'genus', 'species']
    x = np.arange(len(levels))
    width = 0.2

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, method in enumerate(methods):
        if method not in acc_df.columns: continue
        vals = [acc_df[acc_df['level']==lvl][method].values[0]
                if lvl in acc_df['level'].values else 0 for lvl in levels]
        bars = ax.bar(x + i*width, vals, width, label=method, color=get_color(method))
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                    f'{v:.1f}%', ha='center', fontsize=8)

    ax.set_xticks(x + width*1.5)
    ax.set_xticklabels(['Family', 'Genus', 'Species'])
    ax.set_ylabel('Accuracy (%)')
    ax.set_title('Classification Accuracy by Taxonomic Level')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3)
    save_default(fig, outpath)

def plot_error_distribution(error_df, outpath):
    """图11: 错误类型分布堆叠柱状图"""
    set_style()
    methods = error_df['method'].unique()
    error_types = ['correct', 'under_classified', 'mis_classified', 'unassigned']

    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(methods))
    bottom = np.zeros(len(methods))
    colors = {'correct': '#2ECC71', 'under_classified': '#F39C12',
              'mis_classified': '#E74C3C', 'unassigned': '#95A5A6'}

    for etype in error_types:
        if etype not in error_df.columns: continue
        vals = [error_df[error_df['method']==m][etype].values[0]
                if m in error_df['method'].values else 0 for m in methods]
        ax.bar(x, vals, bottom=bottom, label=etype.replace('_', ' ').title(),
               color=colors[etype])
        bottom += np.array(vals)

    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_ylabel('Proportion (%)')
    ax.set_title('Classification Error Type Distribution (Species Level)')
    ax.legend(loc='upper right')
    save_default(fig, outpath)

def plot_coverage_stratified(stratified_df, outpath):
    """图12: 覆盖度梯度 × 分类准确率折线图"""
    set_style()
    fig, ax = plt.subplots(figsize=(9, 6))
    levels = ['family', 'genus', 'species']
    cov_order = ['full_length', 'near_full', 'partial', 'fragment', 'highly_fragmented']
    cov_labels = ['100%', '80%', '60%', '40%', '20%']
    markers = ['o', 's', '^']

    for level, marker in zip(levels, markers):
        ld = stratified_df[stratified_df['level'] == level]
        xs = [i for i, c in enumerate(cov_order) if c in ld['stratum'].values]
        ys = [ld[ld['stratum']==c]['accuracy'].values[0] if c in ld['stratum'].values else np.nan for c in cov_order]
        ax.plot(xs, ys, marker=marker, linewidth=2, label=level.title())

    ax.set_xticks(range(len(cov_labels)))
    ax.set_xticklabels(cov_labels)
    ax.set_xlabel('Genome Coverage')
    ax.set_ylabel('Accuracy (%)')
    ax.set_title('Classification Accuracy by Coverage Level')
    ax.legend()
    ax.grid(alpha=0.3)
    save_default(fig, outpath)

def main():
    parser = argparse.ArgumentParser('评估四绘图')
    parser.add_argument('--accuracy-tsv', help='分类准确率汇总 TSV')
    parser.add_argument('--error-tsv', help='错误类型分布 TSV')
    parser.add_argument('--stratified-tsv', help='分层准确率 TSV')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.accuracy_tsv and os.path.exists(args.accuracy_tsv):
        adf = pd.read_csv(args.accuracy_tsv, sep='\t')
        plot_accuracy_bar(adf, os.path.join(args.outdir, 'Fig10_Classification_Accuracy'))

    if args.error_tsv and os.path.exists(args.error_tsv):
        edf = pd.read_csv(args.error_tsv, sep='\t')
        plot_error_distribution(edf, os.path.join(args.outdir, 'Fig11_Error_Types'))

    if args.stratified_tsv and os.path.exists(args.stratified_tsv):
        sdf = pd.read_csv(args.stratified_tsv, sep='\t')
        plot_coverage_stratified(sdf, os.path.join(args.outdir, 'Fig12_Coverage_vs_Accuracy'))

if __name__ == '__main__':
    main()
