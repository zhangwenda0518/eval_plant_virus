#!/usr/bin/env python3
"""宿主过滤消融实验 — SCI配色 小提琴+箱线+折线三合一 v2"""

import argparse, os
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F']
GROUP_ORDER = ['D0_Baseline', 'D1_Kraken2', 'D2_HISAT2', 'D3_K2+HS2', 'D4_Full']


def plot_violin_box(df, outpath):
    """图A: 小提琴+箱线 病毒保留率（seaborn 自动处理空组）"""
    retention = df[df['Group'] != 'D0_Baseline'][['Group', 'Virus_Retention']].dropna()
    # 过滤全0或常数组
    valid_groups = []
    for g in GROUP_ORDER[1:]:
        vals = retention[retention['Group'] == g]['Virus_Retention']
        if len(vals) > 2 and vals.std() > 0.01:
            valid_groups.append(g)
    retention = retention[retention['Group'].isin(valid_groups)]
    if len(retention) == 0:
        print("  [WARN] No valid violin data, skipping Fig_A")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.violinplot(data=retention, x='Group', y='Virus_Retention',
                   order=valid_groups, palette=PALETTE[1:1+len(valid_groups)],
                   inner='box', cut=0, ax=ax)
    ax.set_ylabel('Virus Read Retention (%)', fontsize=13)
    ax.set_title('Host Depletion — Virus Retention', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    sns.despine()
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_line_enrichment(summary_df, outpath):
    """图B: 折线图 — 去除率 + 保留率 + 富集倍数"""
    fig, ax1 = plt.subplots(figsize=(8, 6))
    x = range(len(summary_df))
    labels = summary_df['Group'].tolist()

    removal = summary_df['Host_Depletion_Rate(%)'].values
    retention = summary_df['Virus_Retention_Rate(%)'].values
    enrichment = summary_df['Virus_Enrichment_Fold'].values

    ax1.plot(x, removal, 'o-', color=PALETTE[0], linewidth=2.5, markersize=9, label='Host Removal (%)')
    ax1.plot(x, retention, 's-', color=PALETTE[1], linewidth=2.5, markersize=9, label='Virus Retention (%)')
    ax1.set_ylabel('Rate (%)', fontsize=13)
    ax1.set_ylim(0, 110)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=9, rotation=10)

    ax2 = ax1.twinx()
    ax2.plot(x, enrichment, 'D-', color=PALETTE[2], linewidth=2.5, markersize=9, label='Enrichment (×)')
    ax2.set_ylabel('Virus Enrichment (fold)', fontsize=13, color=PALETTE[2])
    ax2.tick_params(axis='y', labelcolor=PALETTE[2])

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1+lines2, labels1+labels2, loc='center right', fontsize=10, frameon=True)
    ax1.set_title('Host Depletion Ablation Study', fontsize=14, fontweight='bold')
    ax1.grid(alpha=0.3)
    sns.despine(right=False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_step_waterfall(detail_df, outpath):
    """图C: 瀑布流 — 各阶段 reads 中位数变化"""
    stages = GROUP_ORDER
    med_host = [detail_df[detail_df['Group']==s]['Host_Reads'].median() for s in stages]
    med_virus = [detail_df[detail_df['Group']==s]['Virus_Reads'].median() for s in stages]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(stages))
    ax.bar(x, med_virus, color=PALETTE[2], alpha=0.8, label='Virus Reads')
    ax.bar(x, med_host, bottom=med_virus, color=PALETTE[0], alpha=0.8, label='Host Reads')

    total = [h+v for h, v in zip(med_host, med_virus)]
    for i, (mv, t) in enumerate(zip(med_virus, total)):
        ax.annotate(f'V:{mv:,.0f}', (i, mv/2), ha='center', fontsize=8, color='white', fontweight='bold')
        ax.annotate(f'T:{t:,.0f}', (i, t), textcoords="offset points", xytext=(0, 5),
                    ha='center', fontsize=8, color='#333')

    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=10, rotation=10)
    ax.set_ylabel('Median Reads per Sample', fontsize=13)
    ax.set_title('Step-wise Host Depletion (Median Reads)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3)
    sns.despine()
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def main():
    parser = argparse.ArgumentParser('宿主过滤绘图 v2')
    parser.add_argument('--detail-tsv', required=True)
    parser.add_argument('--summary-tsv', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    detail = pd.read_csv(args.detail_tsv, sep='\t')
    summary = pd.read_csv(args.summary_tsv, sep='\t')

    plot_violin_box(detail, os.path.join(args.outdir, 'Fig_host_A_violin_box.png'))
    plot_line_enrichment(summary, os.path.join(args.outdir, 'Fig_host_B_line.png'))
    plot_step_waterfall(detail, os.path.join(args.outdir, 'Fig_host_C_waterfall.png'))

if __name__ == '__main__':
    main()
