#!/usr/bin/env python3
"""评估一：已知病毒检测方法比较 — 绘图"""

import argparse, os, sys
import pandas as pd, numpy as np
from plot_utils import set_style, get_color, save_default

def plot_f1_vs_lod(eval_df, outpath):
    """图1: F1 × LoD因子 折线图"""
    set_style()
    fig, ax = plt.subplots(figsize=(9, 6))

    for method in eval_df['method'].unique():
        md = eval_df[eval_df['method'] == method]
        grouped = md.groupby('abundance')['f1'].agg(['mean', 'std'])
        ax.errorbar(grouped.index, grouped['mean'], yerr=grouped['std'],
                    label=method, color=get_color(method), marker='o',
                    linewidth=2, markersize=7, capsize=4)

    ax.set_xscale('log')
    ax.set_xlabel('Virus Abundance (LoD factor)')
    ax.set_ylabel('F1 Score')
    ax.set_title('Detection Performance vs Virus Abundance')
    ax.legend(loc='lower right')
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.05)
    save_default(fig, outpath)

def plot_pr_scatter(eval_df, outpath):
    """图2: Precision × Recall 散点图"""
    set_style()
    fig, ax = plt.subplots(figsize=(8, 8))

    for method in eval_df['method'].unique():
        md = eval_df[eval_df['method'] == method]
        avg = md.groupby('abundance')[['precision', 'recall']].mean()
        ax.scatter(avg['recall'], avg['precision'], label=method,
                   color=get_color(method), s=80, alpha=0.8)
        # 连线
        ax.plot(avg['recall'], avg['precision'], color=get_color(method),
                alpha=0.4, linewidth=1)

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Trade-off Across Abundance Levels')
    ax.legend(loc='lower left')
    ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    save_default(fig, outpath)

def plot_runtime(eval_df, outpath):
    """图3: 运行时间对比柱状图"""
    set_style()
    if 'runtime' not in eval_df.columns:
        print("[WARN] No runtime column, skipping runtime plot")
        return
    fig, ax = plt.subplots(figsize=(9, 6))
    methods = eval_df['method'].unique()
    times = [eval_df[eval_df['method']==m]['runtime'].mean() for m in methods]
    colors = [get_color(m) for m in methods]
    bars = ax.bar(methods, times, color=colors, edgecolor='white')
    ax.set_ylabel('Runtime (min/sample)')
    ax.set_title('Average Runtime per Sample')
    for bar, t in zip(bars, times):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.1,
                f'{t:.1f}', ha='center', fontsize=10)
    save_default(fig, outpath)

def main():
    parser = argparse.ArgumentParser('评估一绘图')
    parser.add_argument('--eval-tsv', required=True, help='evaluation_metrics.tsv')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.eval_tsv, sep='\t')
    plot_f1_vs_lod(df, os.path.join(args.outdir, 'Fig1_F1_vs_LoD'))
    plot_pr_scatter(df, os.path.join(args.outdir, 'Fig2_PR_Scatter'))
    plot_runtime(df, os.path.join(args.outdir, 'Fig3_Runtime'))

if __name__ == '__main__':
    main()
