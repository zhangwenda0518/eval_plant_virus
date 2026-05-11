#!/usr/bin/env python3
"""宿主过滤消融实验 — 绘图"""

import argparse, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt
from plot_utils import set_style, get_color, save_default

def plot_depletion_ablation(summary_tsv, outpath):
    """图13: 宿主去除率 & 病毒保留率 双Y轴柱状图"""
    set_style()
    df = pd.read_csv(summary_tsv, sep='\t')
    configs = df['config'].tolist()
    host_removal = df['host_removal_pct'].tolist()
    virus_retention = df['virus_retention_pct'].tolist()
    gain = df['abundance_gain'].tolist()

    fig, ax1 = plt.subplots(figsize=(10, 6))
    x = np.arange(len(configs))
    width = 0.35

    bars1 = ax1.bar(x - width/2, host_removal, width, label='Host Removal (%)',
                     color='#E74C3C', alpha=0.8)
    ax1.set_ylabel('Host Removal (%)', color='#E74C3C')
    ax1.set_ylim(0, 105)

    ax2 = ax1.twinx()
    bars2 = ax2.bar(x + width/2, virus_retention, width, label='Virus Retention (%)',
                     color='#2ECC71', alpha=0.8)
    ax2.set_ylabel('Virus Retention (%)', color='#2ECC71')
    ax2.set_ylim(0, 105)

    ax1.set_xticks(x)
    ax1.set_xticklabels(configs)
    ax1.set_xlabel('Filtering Configuration')
    ax1.set_title('Host Depletion Ablation Study')

    # 丰度增益标注
    for i, (xi, g) in enumerate(zip(x, gain)):
        ax1.annotate(f'×{g:.1f}', (xi, host_removal[i]+1),
                     ha='center', fontsize=9, color='#333333')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1+lines2, labels1+labels2, loc='lower right')
    save_default(fig, outpath)

def main():
    parser = argparse.ArgumentParser('宿主过滤绘图')
    parser.add_argument('--summary-tsv', required=True, help='消融实验汇总 TSV')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plot_depletion_ablation(args.summary_tsv, os.path.join(args.outdir, 'Fig13_Host_Depletion'))

if __name__ == '__main__':
    main()
