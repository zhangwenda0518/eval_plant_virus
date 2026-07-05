#!/usr/bin/env python3
"""
评估四：病毒分类方法比较 — 汇总绘图
从 calc_classification_stratified.py 输出的 stratified_accuracy.tsv 生成对比图。

用法:
  python plot_classification_comparison.py \
      --input step9_classification/analysis/ \
      --outdir step9_classification/analysis/plots/
"""

import argparse, os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

PALETTE = ['#4C72B0', '#55A868', '#C44E52', '#7E6148', '#E64B35', '#F39B7F']
METHOD_DIRS = ['mmseqs', 'vitap', 'acvirus', 'vcontact3', 'integrated']
METHOD_LABELS = {'mmseqs': 'MMseqs2', 'vitap': 'VITAP',
                 'acvirus': 'ACVirus', 'vcontact3': 'vConTACT3',
                 'integrated': 'Integrated'}
ALL_LEVELS = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
CORE_LEVELS = ['family', 'genus', 'species']
HIGH_LEVELS = ['realm', 'kingdom', 'phylum', 'class', 'order']


def _safe_accuracy(sub):
    """返回 accuracy (correct/assigned)，NaN 填 0"""
    vals = sub['accuracy'].values
    return np.nan_to_num(vals, nan=0.0)


def load_all(input_dir):
    """加载所有方法的 stratified_accuracy.tsv"""
    df_all = []
    for d in METHOD_DIRS:
        p = os.path.join(input_dir, d, 'stratified_accuracy.tsv')
        if os.path.exists(p):
            df = pd.read_csv(p, sep='\t')
            df['Method'] = METHOD_LABELS.get(d, d)
            df_all.append(df)
    if not df_all:
        print(f"[ERROR] No stratified_accuracy.tsv found in {input_dir}")
        sys.exit(1)
    return pd.concat(df_all)


def plot_combined(df_all, outdir):
    """一页三图：Family/Genus/Species × 覆盖度折线 — 准确率（correct/assigned）"""
    sns.set_theme(style='whitegrid', font_scale=1.1)
    levels = CORE_LEVELS
    methods = list(METHOD_LABELS.values())

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, level in zip(axes, levels):
        sub = df_all[df_all['level'] == level].copy()
        sub = sub[sub['stratum'] != 'OVERALL']
        sub['cov_num'] = sub['stratum'].str.replace('%', '').astype(int)
        sub = sub.sort_values('cov_num')

        for i, method in enumerate(methods):
            md = sub[sub['Method'] == method]
            if md.empty:
                continue
            ax.plot(md['stratum'], _safe_accuracy(md), marker='o', linewidth=2,
                    color=PALETTE[i], label=method, markersize=6)

        ax.set_title(f'{level.capitalize()} Level', fontsize=14, fontweight='bold')
        ax.set_xlabel('Coverage')
        ax.set_ylabel('Accuracy (%) = correct / assigned')
        ax.set_ylim(0, 105)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    fig.suptitle('Classification Accuracy by Coverage Level', fontsize=16,
                 fontweight='bold', y=1.02)
    plt.tight_layout()
    p = os.path.join(outdir, 'Fig_Classification_Comparison.png')
    fig.savefig(p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {p}")


def plot_higher_levels(df_all, outdir):
    """高层级对比图：Realm→Order × 覆盖度"""
    avail_levels = [l for l in HIGH_LEVELS if l in df_all['level'].values]
    if not avail_levels:
        return
    sns.set_theme(style='whitegrid', font_scale=1.0)
    methods = list(METHOD_LABELS.values())
    n = len(avail_levels)
    fig, axes = plt.subplots(1, n, figsize=(5*n, 5))

    for ax, level in zip(axes if n > 1 else [axes], avail_levels):
        sub = df_all[df_all['level'] == level].copy()
        sub = sub[sub['stratum'] != 'OVERALL']
        if sub.empty:
            continue
        sub['cov_num'] = sub['stratum'].str.replace('%', '').astype(int)
        sub = sub.sort_values('cov_num')
        for i, method in enumerate(methods):
            md = sub[sub['Method'] == method]
            if md.empty:
                continue
            ax.plot(md['stratum'], _safe_accuracy(md), marker='o', linewidth=2,
                    color=PALETTE[i], label=method, markersize=5)
        ax.set_title(f'{level.capitalize()}', fontsize=13, fontweight='bold')
        ax.set_xlabel('Coverage'); ax.set_ylabel('Accuracy (%)')
        ax.set_ylim(0, 105)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        ax.tick_params(axis='x', rotation=45)

    fig.suptitle('High-Level Classification Accuracy', fontsize=15, fontweight='bold', y=1.03)
    plt.tight_layout()
    p = os.path.join(outdir, 'Fig_Classification_Higher_Levels.png')
    fig.savefig(p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {p}")


def plot_overall_bar(df_all, outdir):
    """Overall 柱状图：全部等级 × 准确率 + 比对率"""
    overall = df_all[df_all['stratum'] == 'OVERALL'].copy()
    if overall.empty:
        return
    sns.set_theme(style='whitegrid', font_scale=1.1)
    level_order = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    levels = [l for l in level_order if l in overall['level'].values]
    if not levels:
        levels = CORE_LEVELS
    methods = [m for m in METHOD_LABELS.values()
               if m in overall['Method'].values]
    x = np.arange(len(levels))
    width = 0.15
    n_methods = len(methods)

    fig, ax = plt.subplots(figsize=(12, 6))
    for i, method in enumerate(methods):
        md = overall[overall['Method'] == method]
        vals = []
        for lvl in levels:
            v = md[md['level'] == lvl]['accuracy'].values
            vals.append(v[0] if len(v) > 0 and not pd.isna(v[0]) else 0)
        offset = (i - n_methods/2 + 0.5) * width
        bars = ax.bar(x + offset, vals, width, label=method,
                      color=PALETTE[i % len(PALETTE)])
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.5, f'{v:.1f}', ha='center',
                        fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels([l.capitalize() for l in levels])
    ax.set_ylabel('Accuracy (%) = correct / assigned')
    ax.set_title('Overall Classification Accuracy')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3)
    p = os.path.join(outdir, 'Fig_Classification_Overall_Bar.png')
    fig.savefig(p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {p}")


def plot_unassigned(df_all, outdir):
    """缺失率折线图（使用 missing_rate 列）"""
    sub = df_all[df_all['stratum'] != 'OVERALL'].copy()
    if sub.empty:
        return
    sub['cov_num'] = sub['stratum'].str.replace('%', '').astype(int)
    sub = sub.sort_values('cov_num')
    levels = ['family', 'genus', 'species']

    # 优先使用 missing_rate，否则从 unassigned/total 计算
    if 'missing_rate' in sub.columns:
        sub['miss_val'] = sub['missing_rate']
    elif 'unassigned' in sub.columns:
        sub['miss_val'] = sub['unassigned'] / sub['total'] * 100
    else:
        return

    sns.set_theme(style='whitegrid', font_scale=1.1)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, level in zip(axes, levels):
        ld = sub[sub['level'] == level]
        for i, method in enumerate(METHOD_LABELS.values()):
            md = ld[ld['Method'] == method]
            if md.empty:
                continue
            ax.plot(md['stratum'], md['miss_val'].values, marker='s', linewidth=2,
                    color=PALETTE[i], label=method)
        ax.set_title(f'Missing Rate — {level.capitalize()}')
        ax.set_xlabel('Coverage'); ax.set_ylabel('Missing (%)')
        ax.set_ylim(0, max(25, ax.get_ylim()[1]))
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    plt.tight_layout()
    p = os.path.join(outdir, 'Fig_Missing_Rate.png')
    fig.savefig(p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {p}")


def plot_assignment_rate(df_all, outdir):
    """比对率折线图（使用 assignment_rate 列）"""
    sub = df_all[df_all['stratum'] != 'OVERALL'].copy()
    if 'assignment_rate' not in sub.columns or sub.empty:
        return
    sub['cov_num'] = sub['stratum'].str.replace('%', '').astype(int)
    sub = sub.sort_values('cov_num')
    levels = ['family', 'genus', 'species']

    sns.set_theme(style='whitegrid', font_scale=1.1)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, level in zip(axes, levels):
        ld = sub[sub['level'] == level]
        for i, method in enumerate(METHOD_LABELS.values()):
            md = ld[ld['Method'] == method]
            if md.empty:
                continue
            ax.plot(md['stratum'], md['assignment_rate'].values, marker='D', linewidth=2,
                    color=PALETTE[i], label=method)
        ax.set_title(f'Assignment Rate — {level.capitalize()}')
        ax.set_xlabel('Coverage'); ax.set_ylabel('Assignment (%)')
        ax.set_ylim(0, 105)
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    plt.tight_layout()
    p = os.path.join(outdir, 'Fig_Assignment_Rate.png')
    fig.savefig(p, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {p}")


def main():
    parser = argparse.ArgumentParser(description='分类评估汇总绘图')
    parser.add_argument('--input', default='step9_classification/analysis',
                        help='含 mmseqs/vitap/acvirus/integrated 子目录的 analysis 目录')
    parser.add_argument('--outdir', required=True, help='输出目录')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df_all = load_all(args.input)
    print(f"Loaded {len(df_all)} rows from {df_all['Method'].nunique()} methods")

    # 导出合并 TSV
    merged = os.path.join(args.outdir, 'classification_all_methods.tsv')
    df_all.to_csv(merged, sep='\t', index=False, na_rep='NA')
    print(f"  Merged TSV: {merged}")

    # 图1: 核心三级 (Family/Genus/Species) × 覆盖度 — 准确率
    plot_combined(df_all, args.outdir)

    # 图2: 高层级 (Realm→Order) × 覆盖度 — 准确率
    plot_higher_levels(df_all, args.outdir)

    # 图3: Overall 柱状图 (全部等级) — 准确率
    plot_overall_bar(df_all, args.outdir)

    # 图4: 缺失率
    plot_unassigned(df_all, args.outdir)

    # 图5: 比对率
    plot_assignment_rate(df_all, args.outdir)

    # ── 汇总表格 ──
    level_order = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    avail = [l for l in level_order if l in df_all['level'].values]
    for level in (avail if avail else ['family', 'genus', 'species']):
        sub = df_all[df_all['level'] == level].copy()

        # Table 1: Overall × Method (accuracy, assignment_rate, missing_rate)
        overall_sub = sub[sub['stratum'] == 'OVERALL']
        val_cols = ['accuracy', 'assignment_rate', 'missing_rate', 'correct', 'total', 'assigned', 'unassigned']
        val_cols = [c for c in val_cols if c in overall_sub.columns]
        overall_tbl = overall_sub.pivot_table(
            index='Method', values=val_cols)
        overall_path = os.path.join(args.outdir, f'table_{level}_overall.tsv')
        overall_tbl.to_csv(overall_path, sep='\t', na_rep='NA')
        print(f"  Table: {overall_path}")
        print(f"\n  === {level.upper()} (Overall) ===")
        print(overall_tbl.to_string(float_format=lambda x: f'{x:.1f}'))

        # Table 2: Coverage × Method (accuracy)
        cov_sub = sub[sub['stratum'] != 'OVERALL'].copy()
        if not cov_sub.empty:
            cov_pivot = cov_sub.pivot_table(
                index='stratum', columns='Method', values='accuracy',
                aggfunc='first')
            cov_pivot['_sort'] = cov_pivot.index.str.replace('%', '').astype(int)
            cov_pivot = cov_pivot.sort_values('_sort', ascending=False).drop(columns='_sort')
            cov_path = os.path.join(args.outdir, f'table_{level}_by_coverage.tsv')
            cov_pivot.to_csv(cov_path, sep='\t', na_rep='NA')
            print(f"\n  === {level.upper()} (By Coverage) ===")
            print(cov_pivot.to_string(float_format=lambda x: f'{x:.1f}'))

    print(f"\nDone. All output: {args.outdir}")


if __name__ == '__main__':
    main()
