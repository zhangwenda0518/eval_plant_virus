#!/usr/bin/env python3
"""
节段病毒分类评估 — 比对率、准确率、缺失率 + 工具间一致性

用法:
  python eval_segmented_virus.py \
      --selected step1_eval_viruses/selected_viruses.tsv \
      --meta step4_classification_eval/test_metadata_full.tsv \
      --analysis step9_classification/analysis/ \
      --outdir step9_classification/analysis/segmented/
"""

import argparse, os, sys
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

TAX_LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
PALETTE = ['#4C72B0', '#55A868', '#C44E52', '#7E6148', '#E64B35']
METHODS = {
    'acvirus': 'ACVirus',
    'mmseqs': 'MMseqs2',
    'vitap': 'VITAP',
    'vcontact3': 'vConTACT3',
    'integrated': 'Integrated',
}


def load_data(args):
    """加载所有数据并过滤节段病毒序列"""
    # 节段病毒 accession 列表
    sel = pd.read_csv(args.selected, sep='\t')
    seg_acc = set(sel.loc[sel['genome_type'].str.contains('Segmented', case=False, na=False),
                          'accession'].unique())
    print(f"Segmented virus accessions: {len(seg_acc)}")

    # 元数据
    meta = pd.read_csv(args.meta, sep='\t')
    # 过滤节段病毒
    meta_seg = meta[meta['source_accession'].isin(seg_acc)].copy()
    print(f"Segmented virus test sequences: {len(meta_seg)} / {len(meta)}")

    # 各工具预测
    predictions = {}
    for key, label in METHODS.items():
        path = os.path.join(args.analysis, '..', 'integrated', f'standardized_{key}.tsv')
        if os.path.exists(path):
            predictions[label] = pd.read_csv(path, sep='\t')
    print(f"Loaded predictions for: {list(predictions.keys())}")
    return meta_seg, predictions


def calc_stratified_metrics(pred_df, meta_df, level, pred_map):
    """计算节段病毒的比对率/准确率/缺失率（按覆盖度分层）"""
    results = defaultdict(lambda: {"correct": 0, "total": 0, "unassigned": 0})

    for _, row in meta_df.iterrows():
        sid = row.get('seq_id') or row.get('source_accession') or row.get('id')
        cov = int(row.get('coverage_pct', 100))
        stratum = f"{cov}%"
        true_val = str(row.get(level, row.get(level.lower(), ''))).strip()

        pred_val = pred_map.get(sid, None)

        results[stratum]["total"] += 1
        if pred_val is None or pred_val in ('', 'Unassigned', 'none', 'NA', 'nan'):
            results[stratum]["unassigned"] += 1
        elif pred_val == true_val:
            results[stratum]["correct"] += 1

    rows = []
    strata_order = sorted(results.keys(), key=lambda x: -int(x.replace('%', '')))
    for stratum in strata_order:
        r = results[stratum]
        assigned = r['total'] - r['unassigned']
        acc = r['correct'] / assigned * 100 if assigned > 0 else None
        assign_rate = assigned / r['total'] * 100 if r['total'] > 0 else 0
        miss_rate = r['unassigned'] / r['total'] * 100 if r['total'] > 0 else 0
        rows.append({
            'level': level, 'stratum': stratum,
            'total': r['total'], 'assigned': assigned,
            'unassigned': r['unassigned'], 'correct': r['correct'],
            'accuracy': round(acc, 2) if acc is not None else None,
            'assignment_rate': round(assign_rate, 2),
            'missing_rate': round(miss_rate, 2),
        })

    # Overall
    total_c = sum(r['correct'] for r in results.values())
    total_a = sum(r['total'] for r in results.values())
    total_u = sum(r['unassigned'] for r in results.values())
    assigned = total_a - total_u
    overall_acc = total_c / assigned * 100 if assigned > 0 else None
    overall_assign = assigned / total_a * 100 if total_a > 0 else 0
    overall_miss = total_u / total_a * 100 if total_a > 0 else 0
    rows.append({
        'level': level, 'stratum': 'OVERALL',
        'total': total_a, 'assigned': assigned,
        'unassigned': total_u, 'correct': total_c,
        'accuracy': round(overall_acc, 2) if overall_acc is not None else None,
        'assignment_rate': round(overall_assign, 2),
        'missing_rate': round(overall_miss, 2),
    })
    return rows


def compute_per_tool_metrics(meta_seg, predictions, outdir):
    """每个工具分别计算分层指标"""
    all_rows = []
    for method, pred_df in predictions.items():
        # 构建 pred_map
        pred_map = {}
        for _, row in pred_df.iterrows():
            sid = row.get('contig_id') or row.get('seq_id') or row.get('id')
            if not sid:
                continue
            for lvl in TAX_LEVELS:
                val = row.get(lvl, row.get(lvl.lower(), ''))
                if val and str(val).strip() not in ('', 'NA', 'none', 'Unassigned', 'nan'):
                    if sid not in pred_map:
                        pred_map[sid] = {}
                    pred_map[sid][lvl] = str(val).strip()

        for level in TAX_LEVELS:
            # 检查 metadata 中是否有该等级的 ground truth
            has_truth = level in meta_seg.columns or level.lower() in meta_seg.columns
            if not has_truth:
                continue
            lvl_map = {sid: d[level] for sid, d in pred_map.items() if level in d}
            rows = calc_stratified_metrics(pred_df, meta_seg, level, lvl_map)
            for r in rows:
                r['method'] = method
            all_rows.extend(rows)

    out_df = pd.DataFrame(all_rows)
    path = os.path.join(outdir, 'segmented_stratified_accuracy.tsv')
    out_df.to_csv(path, sep='\t', index=False, na_rep='NA')
    print(f"Saved: {path}")
    return out_df


def compute_inter_tool_consistency(meta_seg, predictions, outdir):
    """计算工具间一致性：对于同一个序列，多个工具的分类是否一致"""
    # 构建 序列 × 工具 × 等级 的矩阵
    # seq_tool_level[seq_id][tool][level] = classification_value
    seq_data = defaultdict(lambda: defaultdict(dict))

    for method, pred_df in predictions.items():
        for _, row in pred_df.iterrows():
            sid = row.get('contig_id') or row.get('seq_id') or row.get('id')
            if not sid:
                continue
            for lvl in TAX_LEVELS:
                val = row.get(lvl, row.get(lvl.lower(), ''))
                if val and str(val).strip() not in ('', 'NA', 'none', 'Unassigned', 'nan'):
                    seq_data[sid][method][lvl] = str(val).strip()

    # 过滤到节段病毒序列
    seg_ids = set()
    for _, row in meta_seg.iterrows():
        sid = row.get('seq_id') or row.get('source_accession') or row.get('id')
        seg_ids.add(sid)

    consistency_rows = []
    for level in TAX_LEVELS:
        # 只分析 ≥2 个工具都分类了的序列
        # 注意：不包括 Integrated（因为它本身就是合成的）
        real_tools = [m for m in predictions.keys() if m != 'Integrated']
        agree_count = 0
        disagree_count = 0
        partial_count = 0  # 只有1个工具分类了
        none_count = 0     # 没有任何工具分类

        for sid in seg_ids:
            if sid not in seq_data:
                none_count += 1
                continue

            vals = {}
            for tool in real_tools:
                if tool in seq_data[sid] and level in seq_data[sid][tool]:
                    vals[tool] = seq_data[sid][tool][level]

            n_tools = len(vals)
            if n_tools == 0:
                none_count += 1
            elif n_tools == 1:
                partial_count += 1
            else:
                # 检查所有工具的取值是否一致
                unique_vals = set(vals.values())
                if len(unique_vals) == 1:
                    agree_count += 1
                else:
                    disagree_count += 1

        total = agree_count + disagree_count + partial_count + none_count
        multi_tool = agree_count + disagree_count  # ≥2 工具参与
        consistency_rows.append({
            'level': level,
            'total_sequences': total,
            'all_agree': agree_count,
            'disagree': disagree_count,
            'single_tool_only': partial_count,
            'none_classified': none_count,
            'consistency_rate': round(agree_count / multi_tool * 100, 1) if multi_tool > 0 else None,
            'multi_tool_rate': round(multi_tool / total * 100, 1) if total > 0 else 0,
        })

    out_df = pd.DataFrame(consistency_rows)
    path = os.path.join(outdir, 'segmented_consistency.tsv')
    out_df.to_csv(path, sep='\t', index=False, na_rep='NA')
    print(f"Saved: {path}")
    return out_df


def print_summary(metrics_df, consistency_df, meta_seg):
    """打印摘要"""
    print(f"\n{'='*70}")
    print(f"  节段病毒分类评估摘要")
    print(f"{'='*70}")
    print(f"  节段病毒测试序列: {len(meta_seg)}")
    print()

    # 核心等级 Overall 结果
    core_levels = ['Family', 'Genus', 'Species']
    overall = metrics_df[metrics_df['stratum'] == 'OVERALL']
    for level in core_levels:
        sub = overall[overall['level'] == level]
        if sub.empty:
            continue
        print(f"--- {level.upper()} (Overall) ---")
        print(f"  {'Method':12s} | {'Assign':>8s} | {'Accuracy':>8s} | {'Missing':>8s}")
        print(f"  {'-'*12}-+{'-'*10}+{'-'*10}+{'-'*10}")
        for _, r in sub.sort_values('method').iterrows():
            acc_str = f"{r['accuracy']:.1f}%" if not pd.isna(r['accuracy']) else "N/A"
            print(f"  {r['method']:12s} | {r['assignment_rate']:5.1f}%   | {acc_str:>7s} | {r['missing_rate']:5.1f}%")
        print()

    # 一致性
    print(f"--- 工具间一致性 ---")
    print(f"  {'Level':10s} | {'AllAgree':>9s} | {'Disagree':>9s} | {'Single':>7s} | {'None':>5s} | {'Consist':>7s}")
    print(f"  {'-'*10}-+-{'-'*9}-+-{'-'*9}-+-{'-'*7}-+-{'-'*5}-+-{'-'*7}")
    for _, r in consistency_df.iterrows():
        consist_str = f"{r['consistency_rate']:.1f}%" if not pd.isna(r['consistency_rate']) else "N/A"
        print(f"  {r['level']:10s} | {int(r['all_agree']):5d}      | {int(r['disagree']):5d}      | {int(r['single_tool_only']):5d} | {int(r['none_classified']):4d} | {consist_str:>6s}")
    print()


def plot_segmented_results(metrics_df, consistency_df, meta_seg, outdir):
    """节段病毒分类评估综合图"""

    methods_order = list(METHODS.values())
    available_methods = [m for m in methods_order if m in metrics_df['method'].values]
    available_levels = TAX_LEVELS
    metrics_levels_lower = set(str(l).lower() for l in metrics_df['level'].unique())
    available_levels = [l for l in TAX_LEVELS if l.lower() in metrics_levels_lower]

    # ---- 图1: Overall 三指标对比柱状图 (3×1 subplot) ----
    overall = metrics_df[metrics_df['stratum'] == 'OVERALL'].copy()
    if not overall.empty:
        sns.set_theme(style='whitegrid', font_scale=1.05)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        metric_configs = [
            ('accuracy', 'Accuracy (%)', 'Classification Accuracy (correct / assigned)'),
            ('assignment_rate', 'Assignment Rate (%)', 'Assignment Rate (assigned / total)'),
            ('missing_rate', 'Missing Rate (%)', 'Missing Rate (unassigned / total)'),
        ]

        for ax, (col, ylabel, title) in zip(axes, metric_configs):
            sub = overall[overall['level'].isin(available_levels)]

            x = np.arange(len(available_levels))
            n_methods = len(available_methods)
            width = 0.8 / n_methods

            for i, method in enumerate(available_methods):
                md = sub[sub['method'] == method]
                vals = []
                for lvl in available_levels:
                    v = md[md['level'] == lvl][col].values
                    vals.append(v[0] if len(v) > 0 and not pd.isna(v[0]) else 0)
                bars = ax.bar(x + i * width, vals, width, label=method,
                              color=PALETTE[i % len(PALETTE)])
                for bar, v in zip(bars, vals):
                    if v > 0:
                        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                                f'{v:.0f}', ha='center', fontsize=7)

            ax.set_xticks(x + width * (n_methods - 1) / 2)
            ax.set_xticklabels([l.capitalize() for l in available_levels], rotation=45 if len(available_levels) > 5 else 0)
            ax.set_ylabel(ylabel)
            ax.set_title(title, fontsize=12, fontweight='bold')
            ax.set_ylim(0, 110)
            ax.legend(fontsize=7, loc='lower right')
            ax.grid(axis='y', alpha=0.3)

        fig.suptitle('Segmented Virus Classification — Overall Metrics',
                     fontsize=14, fontweight='bold', y=1.01)
        plt.tight_layout()
        path = os.path.join(outdir, 'Fig_Segmented_Overall_Metrics.png')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {path}")

    # ---- 图2: 工具间一致性堆叠柱状图 ----
    if not consistency_df.empty:
        sns.set_theme(style='whitegrid', font_scale=1.1)
        cons = consistency_df.copy()
        cons['level'] = pd.Categorical(cons['level'],
                                       categories=TAX_LEVELS,
                                       ordered=True)
        cons = cons.sort_values('level')

        fig, ax = plt.subplots(figsize=(12, 6))
        cons.set_index('level')[['all_agree', 'disagree', 'single_tool_only', 'none_classified']].plot(
            kind='bar', stacked=True, ax=ax,
            color=['#55A868', '#C44E52', '#F39B7F', '#CCCCCC'])
        ax.set_title('Inter-Tool Consistency — Segmented Viruses', fontsize=14, fontweight='bold')
        ax.set_xlabel('Taxonomic Level')
        ax.set_ylabel('Number of Sequences')
        ax.legend(['All Agree', 'Disagree', 'Single Tool', 'None'], fontsize=9,
                  loc='upper right')
        ax.set_xticklabels([str(l).capitalize() for l in cons['level']], rotation=0)
        ax.grid(axis='y', alpha=0.3)

        # 在 bar 顶部标注一致率
        for i, (_, r) in enumerate(cons.iterrows()):
            if r['consistency_rate'] == r['consistency_rate']:  # not NaN
                total_h = r['all_agree'] + r['disagree'] + r['single_tool_only'] + r['none_classified']
                ax.text(i, total_h + 0.5, f"{r['consistency_rate']:.0f}%",
                        ha='center', fontsize=8, fontweight='bold', color='#55A868')

        plt.tight_layout()
        path = os.path.join(outdir, 'Fig_Segmented_Consistency.png')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {path}")

    # ---- 图3: 核心三级覆盖度折线图 (Accuracy) ----
    core_levels = ['Family', 'Genus', 'Species']
    core_sub = metrics_df[
        (metrics_df['stratum'] != 'OVERALL') &
        (metrics_df['level'].isin(core_levels))
    ].copy()
    if not core_sub.empty:
        core_sub['cov_num'] = core_sub['stratum'].str.replace('%', '').astype(int)
        core_sub = core_sub.sort_values('cov_num')

        sns.set_theme(style='whitegrid', font_scale=1.1)
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        for ax, level in zip(axes, core_levels):
            ld = core_sub[core_sub['level'] == level]
            for i, method in enumerate(available_methods):
                md = ld[ld['method'] == method]
                if md.empty:
                    continue
                vals = np.nan_to_num(md['accuracy'].values, nan=0.0)
                ax.plot(md['stratum'], vals, marker='o', linewidth=2,
                        color=PALETTE[i % len(PALETTE)], label=method, markersize=5)
            ax.set_title(f'{level.capitalize()} Level', fontsize=13, fontweight='bold')
            ax.set_xlabel('Coverage')
            ax.set_ylabel('Accuracy (%)')
            ax.set_ylim(0, 105)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)

        fig.suptitle('Segmented Virus — Accuracy by Coverage Level',
                     fontsize=15, fontweight='bold', y=1.02)
        plt.tight_layout()
        path = os.path.join(outdir, 'Fig_Segmented_Accuracy_by_Coverage.png')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {path}")

    print(f"All plots saved to: {outdir}")


def main():
    parser = argparse.ArgumentParser(description='节段病毒分类评估')
    parser.add_argument('--selected', required=True, help='selected_viruses.tsv')
    parser.add_argument('--meta', required=True, help='test_metadata_full.tsv')
    parser.add_argument('--analysis', default='step9_classification/analysis/',
                        help='analysis 目录 (含各工具 stratified_accuracy.tsv)')
    parser.add_argument('--outdir', required=True, help='输出目录')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta_seg, predictions = load_data(args)

    if len(meta_seg) == 0:
        print("No segmented virus sequences found in metadata!")
        sys.exit(1)

    # 1. 各工具分级指标
    metrics_df = compute_per_tool_metrics(meta_seg, predictions, args.outdir)

    # 2. 工具间一致性
    consistency_df = compute_inter_tool_consistency(meta_seg, predictions, args.outdir)

    # 3. 打印摘要
    print_summary(metrics_df, consistency_df, meta_seg)

    # 4. 绘图
    plot_segmented_results(metrics_df, consistency_df, meta_seg, args.outdir)

    print(f"Done. All output: {args.outdir}")


if __name__ == '__main__':
    main()
