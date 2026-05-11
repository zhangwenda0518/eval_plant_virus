#!/usr/bin/env python3
"""评估二：病毒组装方法比较 — 绘图"""

import argparse, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from plot_utils import set_style, get_color, save_default

def plot_quality_boxplot(quality_df, outpath):
    """图4: NGA50 / Genome fraction 箱线图"""
    set_style()
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    metrics = [
        ('genome_fraction', 'Genome Fraction (%)', axes[0,0]),
        ('NGA50', 'NGA50 (bp)', axes[0,1]),
        ('misassemblies_per_mbp', 'Misassemblies / Mbp', axes[1,0]),
        ('mismatches_per_100kbp', 'Mismatches / 100 kbp', axes[1,1]),
    ]

    for col, ylabel, ax in metrics:
        if col not in quality_df.columns: continue
        methods = quality_df['tool'].unique()
        plot_data = [quality_df[quality_df['tool']==m][col].dropna().values for m in methods]
        bp = ax.boxplot(plot_data, labels=methods, patch_artist=True)
        for patch, method in zip(bp['boxes'], methods):
            patch.set_facecolor(get_color(method))
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)

    fig.suptitle('Assembly Quality Comparison', fontsize=16, fontweight='bold')
    fig.tight_layout()
    save_default(fig, outpath)

def plot_chimeric_rate(chimeric_summary, outpath):
    """图5: 嵌合率对比柱状图"""
    set_style()
    fig, ax = plt.subplots(figsize=(8, 6))
    methods = list(chimeric_summary.keys())
    rates = list(chimeric_summary.values())
    colors = [get_color(m) for m in methods]
    bars = ax.bar(methods, rates, color=colors, edgecolor='white')
    ax.set_ylabel('Chimeric Rate (%)')
    ax.set_title('Chimeric Contig Rate by Assembler')
    for bar, r in zip(bars, rates):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
                f'{r:.2f}%', ha='center')
    save_default(fig, outpath)

def plot_resource_scatter(resource_df, outpath):
    """图6: 运行时间 vs 内存散点图"""
    set_style()
    fig, ax = plt.subplots(figsize=(8, 7))
    for tool in resource_df['tool'].unique():
        td = resource_df[resource_df['tool'] == tool]
        ax.scatter(td['time_seconds'], td['memory_gb'], label=tool,
                   color=get_color(tool), s=100, alpha=0.8)
    ax.set_xlabel('Runtime (seconds)')
    ax.set_ylabel('Peak Memory (GB)')
    ax.set_title('Resource Consumption Trade-off')
    ax.legend()
    ax.grid(alpha=0.3)
    save_default(fig, outpath)

def main():
    parser = argparse.ArgumentParser('评估二绘图')
    parser.add_argument('--quality-csv', help='benchmark_quality_summary.csv')
    parser.add_argument('--chimeric-tsv', help='chimeric_rate_summary.tsv')
    parser.add_argument('--resource-tsv', help='resource_summary.tsv')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.quality_csv and os.path.exists(args.quality_csv):
        qdf = pd.read_csv(args.quality_csv)
        plot_quality_boxplot(qdf, os.path.join(args.outdir, 'Fig4_Assembly_Quality'))

    if args.chimeric_tsv and os.path.exists(args.chimeric_tsv):
        cdf = pd.read_csv(args.chimeric_tsv, sep='\t')
        cdict = dict(zip(cdf.iloc[:,0], cdf.iloc[:,1]))
        plot_chimeric_rate(cdict, os.path.join(args.outdir, 'Fig5_Chimeric_Rate'))

    if args.resource_tsv and os.path.exists(args.resource_tsv):
        rdf = pd.read_csv(args.resource_tsv, sep='\t')
        plot_resource_scatter(rdf, os.path.join(args.outdir, 'Fig6_Resource'))

if __name__ == '__main__':
    main()
