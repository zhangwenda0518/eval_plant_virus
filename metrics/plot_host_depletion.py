#!/usr/bin/env python3
"""宿主过滤消融实验 — SCI配色 小提琴+箱线+折线三合一"""

import argparse, os
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# SCI 论文级配色 (Nature Reviews 风格)
PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F']
GROUP_ORDER = ['D0_Baseline', 'D1_Kraken2', 'D2_HISAT2', 'D3_K2+HS2', 'D4_Full']


def plot_violin_box(df, outpath):
    """图A: 小提琴+箱线 病毒保留率"""
    retention = df[df['Group'] != 'D0_Baseline'][['Group', 'Virus_Retention']].dropna()

    fig, ax = plt.subplots(figsize=(8, 6))
    parts = ax.violinplot(
        [retention[retention['Group'] == g]['Virus_Retention'].values for g in GROUP_ORDER[1:]],
        positions=range(len(GROUP_ORDER[1:])), showmeans=False, showmedians=False
    )
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(PALETTE[i + 1])
        pc.set_alpha(0.7)

    bp = ax.boxplot(
        [retention[retention['Group'] == g]['Virus_Retention'].values for g in GROUP_ORDER[1:]],
        positions=range(len(GROUP_ORDER[1:])), widths=0.15,
        patch_artist=True, showfliers=False
    )
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor('white')
        patch.set_edgecolor(PALETTE[i + 1])
        patch.set_linewidth(1.5)

    ax.set_xticks(range(len(GROUP_ORDER[1:])))
    ax.set_xticklabels(['Kraken2', 'HISAT2', 'K2+HISAT2', 'Full (K2+HS2+rRNA)'], fontsize=11)
    ax.set_ylabel('Virus Read Retention (%)', fontsize=13)
    ax.set_title('Host Depletion — Virus Retention', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    sns.despine()
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_line_enrichment(summary_df, detail_df, outpath):
    """图B: 折线图 — 宿主去除率 + 病毒保留率 + 富集倍数"""
    fig, ax1 = plt.subplots(figsize=(8, 6))
    x = range(len(GROUP_ORDER))

    removal = summary_df['Host_Removal(%)'].values
    retention = summary_df['Virus_Retention(%)'].values
    enrichment = summary_df['Enrichment(x)'].values

    ax1.plot(x, removal, 'o-', color=PALETTE[0], linewidth=2.5, markersize=9, label='Host Removal (%)')
    ax1.plot(x, retention, 's-', color=PALETTE[1], linewidth=2.5, markersize=9, label='Virus Retention (%)')
    ax1.set_ylabel('Rate (%)', fontsize=13)
    ax1.set_ylim(0, 110)
    ax1.set_xticks(x)
    ax1.set_xticklabels(['D0\nBaseline', 'D1\nKraken2', 'D2\nHISAT2', 'D3\nK2+HS2', 'D4\nFull(K2+HS2+rRNA)'], fontsize=10)

    ax2 = ax1.twinx()
    ax2.plot(x, enrichment, 'D-', color=PALETTE[2], linewidth=2.5, markersize=9, label='Enrichment (×)')
    ax2.set_ylabel('Virus Enrichment (fold)', fontsize=13, color=PALETTE[2])
    ax2.tick_params(axis='y', labelcolor=PALETTE[2])

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right', fontsize=11, frameon=True)

    ax1.set_title('Host Depletion Ablation Study', fontsize=14, fontweight='bold')
    ax1.grid(alpha=0.3)
    sns.despine(right=False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_step_waterfall(detail_df, outpath):
    """图C: 瀑布流 — 各样本在不同过滤阶段的 reads 变化"""
    # 取中位数汇总每个阶段
    stages = ['D0_Baseline', 'D1_Kraken2', 'D2_HISAT2', 'D3_K2+HS2', 'D4_Full']
    med_host = [detail_df[detail_df['Group']==s]['Host_Reads'].median() for s in stages]
    med_virus = [detail_df[detail_df['Group']==s]['Virus_Reads'].median() for s in stages]

    fig, ax = plt.subplots(figsize=(8, 6))
    x = range(len(stages))
    ax.fill_between(x, 0, med_virus, color=PALETTE[2], alpha=0.6, label='Virus Reads')
    ax.fill_between(x, med_virus, [mh + mv for mh, mv in zip(med_host, med_virus)],
                    color=PALETTE[0], alpha=0.6, label='Host Reads')

    total = [h + v for h, v in zip(med_host, med_virus)]
    for i, t in enumerate(total):
        ax.annotate(f'{t:,.0f}', (i, t), textcoords="offset points", xytext=(0, 8),
                    ha='center', fontsize=9, color='#333333')

    ax.set_xticks(x)
    ax.set_xticklabels(['Baseline', 'Kraken2', 'HISAT2', 'K2+HS2', 'Full'], fontsize=11)
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
    parser = argparse.ArgumentParser('宿主过滤绘图')
    parser.add_argument('--detail-tsv', required=True, help='host_depletion_detail.tsv')
    parser.add_argument('--summary-tsv', help='host_depletion_report.tsv')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    detail = pd.read_csv(args.detail_tsv, sep='\t')
    summary = pd.read_csv(args.summary_tsv, sep='\t') if args.summary_tsv and os.path.exists(args.summary_tsv) else None

    if summary is None:
        # 从 detail 计算 summary
        summ = detail.groupby('Group').agg(Host_Total=('Host_Reads', 'sum'), Virus_Total=('Virus_Reads', 'sum')).reset_index()
        bh, bv = summ.iloc[0]['Host_Total'], summ.iloc[0]['Virus_Total']
        summ['Host_Removal(%)'] = ((1 - summ['Host_Total'] / bh) * 100).round(2)
        summ['Virus_Retention(%)'] = ((summ['Virus_Total'] / bv) * 100).round(2)
        br = bv / (bh + bv) if (bh + bv) > 0 else 1
        summ['Enrichment(x)'] = summ.apply(
            lambda r: round((r['Virus_Total'] / (r['Host_Total'] + r['Virus_Total'])) / br, 1)
            if (r['Host_Total'] + r['Virus_Total']) > 0 else 0, axis=1)
        summary = summ

    plot_violin_box(detail, os.path.join(args.outdir, 'Fig_host_A_violin_box.png'))
    plot_line_enrichment(summary, detail, os.path.join(args.outdir, 'Fig_host_B_line.png'))
    plot_step_waterfall(detail, os.path.join(args.outdir, 'Fig_host_C_waterfall.png'))


if __name__ == '__main__':
    main()
