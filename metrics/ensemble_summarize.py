#!/usr/bin/env python3
"""
分类集成优化结果总结 + 可视化（三指标输出 + 4 张对比图）

用法:
  python ensemble_summarize.py \
      --opt-tsv ensemble_optimization.tsv \
      --best-tsv ensemble_best_per_level.tsv \
      --outdir ensemble_opt/
"""

import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

TAX_LEVELS = ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species']
METRIC_KEYS = ['accuracy','assignment_rate','f1_score','specificity','combined_score']
METRIC_LABELS = ['Accuracy','Assign Rate','F1 Score','Specificity','Combined']
SORT_METRICS = [('combined_score','综合得分(F1×抗FP)'),('f1_score','F1-score'),('accuracy','准确率')]
SINGLE_TOOLS = ['ACVirus','CAT','MMseqs2','vConTACT3','VITAP','diamond_lca','genomad','metabuli']

# ══════════════════════════════════════
# 辅助
# ══════════════════════════════════════

def _get_all_tools(opt_df):
    return sorted(opt_df[opt_df['strategy'] == 'single']['subset'].unique())


# ══════════════════════════════════════
# Part 1: 三指标文本输出
# ══════════════════════════════════════

def summarize_three_metrics(opt_df, outdir):
    all_best = []
    for sc, sc_name in SORT_METRICS:
        print(f"\n{'='*70}")
        print(f"  最优组合 — {sc_name}")
        for vt_label, vt_title in [('known','已知病毒'),('novel','新病毒')]:
            vdf = opt_df[opt_df['virus_type'] == vt_label]
            print(f"\n  ── {vt_title} ──")
            print(f"  {'Level':10s} | {'Best':28s} | {'Strategy':12s} | {'Acc':>5s} | {'Assign':>6s} | {'F1':>4s} | {'FP%':>5s} | {'CS':>4s} | N")
            print(f"  {'-'*10}-+-{'-'*28}-+-{'-'*12}-+-{'-'*5}-+-{'-'*6}-+-{'-'*4}-+-{'-'*5}-+-{'-'*4}-+--")
            for lvl in TAX_LEVELS:
                ldf = vdf[vdf['level'] == lvl].sort_values(sc, ascending=False, na_position='last')
                if ldf.empty: continue
                b = ldf.iloc[0].to_dict()
                b['sort_metric'] = sc; all_best.append(b)
                a_s  = f"{b['accuracy']:.0f}%" if pd.notna(b['accuracy']) else "NA"
                f_s  = f"{b['f1_score']:.0f}" if pd.notna(b['f1_score']) else "NA"
                fp_s = f"{b['fp_rate']:.0f}%" if pd.notna(b['fp_rate']) else "NA"
                cs_s = f"{b['combined_score']:.0f}" if pd.notna(b['combined_score']) else "NA"
                ar_s = f"{b['assignment_rate']:.1f}%" if pd.notna(b['assignment_rate']) else "NA"
                print(f"  {lvl:10s} | {b['subset']:28s} | {b['strategy']:12s} | {a_s:>5} | {ar_s:>6} | {f_s:>4} | {fp_s:>5} | {cs_s:>4} | {b['n_tools']:>2d}")
    if all_best:
        pd.DataFrame(all_best).to_csv(os.path.join(outdir,'ensemble_best_per_level.tsv'), sep='\t', index=False, na_rep='NA')


# ══════════════════════════════════════
# Part 2: 可视化
# ══════════════════════════════════════

def plot_compare_bars(opt_df, outdir):
    """单工具三指标对比条形图"""
    singles = opt_df[opt_df['strategy'] == 'single'].copy()
    if singles.empty: return
    tools = _get_all_tools(opt_df)
    metrics_plot = [('accuracy','Accuracy (%)'),('f1_score','F1-score'),('combined_score','Combined Score')]
    for vt_label, vt_title in [('known','Known'),('novel','Novel')]:
        vdf = singles[singles['virus_type'] == vt_label]
        if vdf.empty: continue
        sns.set_theme(style='whitegrid', font_scale=1.0)
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        for ax, (col, ylabel) in zip(axes, metrics_plot):
            x = np.arange(len(TAX_LEVELS))
            width = 0.8 / max(len(tools), 1)
            for i, tool in enumerate(tools):
                vals = []
                for lvl in TAX_LEVELS:
                    r = vdf[(vdf['subset']==tool) & (vdf['level']==lvl)]
                    vals.append(r[col].values[0] if len(r)>0 and pd.notna(r[col].values[0]) else 0)
                ax.bar(x + i*width, vals, width, label=tool[:12], color=PALETTE[i % len(PALETTE)])
            ax.set_xticks(x + width*(len(tools)-1)/2)
            ax.set_xticklabels([l[:3] for l in TAX_LEVELS], fontsize=8)
            ax.set_ylabel(ylabel); ax.set_ylim(0, 105)
            ax.set_title(f'{vt_title} — {ylabel}', fontsize=11, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
        # 图例放标题下方
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=min(6, len(handles)),
                   fontsize=8, bbox_to_anchor=(0.5, 0.97), frameon=False)
        fig.suptitle(f'Tool Comparison ({vt_title})', fontsize=14, fontweight='bold', y=1.05)
        plt.tight_layout()
        path = os.path.join(outdir, f'Fig_Compare_{vt_label}.png')
        fig.savefig(path, dpi=200, bbox_inches='tight'); plt.close()
        print(f"  [plot] {path}")


def plot_performance_curves(opt_df, best_df, outdir):
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    all_tools = _get_all_tools(opt_df)
    for idx, (vt_label, vt_title) in enumerate([('known','Known'),('novel','Novel')]):
        ax = axes[idx]
        for tool in all_tools:
            sub = opt_df[(opt_df['subset']==tool) & (opt_df['strategy']=='single') & (opt_df['virus_type']==vt_label)]
            if sub.empty: continue
            sub = sub.set_index('level').reindex(TAX_LEVELS).reset_index()
            ax.plot(sub['level'], sub['combined_score'], marker='o', linewidth=2, label=tool, alpha=0.7)
        best_sub = best_df[best_df['virus_type']==vt_label] if 'virus_type' in best_df.columns else best_df
        if not best_sub.empty:
            best_sub = best_sub.set_index('level').reindex(TAX_LEVELS).reset_index()
            ax.plot(best_sub['level'], best_sub['combined_score'],
                    marker='*', markersize=12, color='red', linestyle='--', linewidth=3.5, label='Best Ens.')
            for _, row in best_sub.iterrows():
                if pd.notna(row.get('combined_score')):
                    ax.annotate(str(row.get('strategy','')), (row['level'], row['combined_score']),
                                textcoords="offset points", xytext=(0,10), ha='center', fontsize=9, color='darkred')
        ax.set_title(f"Classification - {vt_title}", fontsize=14, weight='bold')
        ax.set_xlabel("Levels"); ax.tick_params(axis='x', rotation=30)
        if idx==0: ax.set_ylabel("Combined Score (F1×(1-FPR))")
        ax.set_ylim(-5,105); ax.legend(fontsize=8, loc='lower left')
    # 公共图例放标题下
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=min(6, len(handles)),
               fontsize=8, bbox_to_anchor=(0.5, 0.99), frameon=False)
    plt.tight_layout()
    path = os.path.join(outdir, 'performance_lift_curves.png')
    plt.savefig(path, dpi=300); plt.close()
    print(f"  [plot] {path}")


def plot_decoy_fpr(opt_df, outdir):
    sns.set_theme(style="white")
    all_tools = _get_all_tools(opt_df)
    sub = opt_df[(opt_df['subset'].isin(all_tools)) & (opt_df['strategy']=='single') &
                 (opt_df['level']=='Species') & (opt_df['virus_type']=='known')]
    if sub.empty: return
    sub = sub.drop_duplicates(subset=['subset']).sort_values('fp_rate', ascending=False)
    plt.figure(figsize=(10, 6))
    colors = sns.color_palette("Reds_r", len(sub))
    bars = plt.bar(sub['subset'], sub['fp_rate'], color=colors, edgecolor='grey', width=0.6)
    for bar in bars:
        h = bar.get_height()
        plt.text(bar.get_x()+bar.get_width()/2.0, h+1, f"{h:.1f}%", ha='center', va='bottom', fontsize=10)
    plt.title("False Positive Rate on Non-Viral Decoy Sequences", fontsize=13, weight='bold')
    plt.xlabel("Classifiers"); plt.ylabel("False Positive Rate (%)")
    plt.ylim(0, max(sub['fp_rate'])*1.15); sns.despine()
    path = os.path.join(outdir, 'decoy_false_positive_rates.png')
    plt.tight_layout(); plt.savefig(path, dpi=300); plt.close()
    print(f"  [plot] {path}")


def get_metrics_row(df, subset, strategy, level, vt):
    sub = df[(df['subset']==subset) & (df['strategy']==strategy) &
             (df['level']==level) & (df['virus_type']==vt)]
    if sub.empty: return [0.0]*len(METRIC_KEYS)
    row = sub.iloc[0]
    fpr = float(row['fp_rate']) if pd.notna(row['fp_rate']) else 0.0
    return [
        float(row['accuracy']) if pd.notna(row['accuracy']) else 0.0,
        float(row['assignment_rate']) if pd.notna(row['assignment_rate']) else 0.0,
        float(row['f1_score']) if pd.notna(row['f1_score']) else 0.0,
        max(0.0, 100.0-fpr),
        float(row['combined_score']) if pd.notna(row['combined_score']) else 0.0,
    ]


TOOL_COLORS = ['#1f77b4','#2ca02c','#ff7f0e','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#d62728']


def draw_radar_grid(opt_df, best_df, vt_label, vt_title, out_path):
    """2×4 雷达图矩阵：每等级一个子图，五边形=5指标，各线=各工具"""
    sns.set_theme(style="white")
    fig, axes = plt.subplots(2, 4, figsize=(22, 12), subplot_kw=dict(projection='polar'))
    axes = axes.flatten()
    num_vars = len(METRIC_KEYS)
    angles = [n / float(num_vars) * 2 * np.pi for n in range(num_vars)]
    angles_pad = angles + angles[:1]

    # 收集所有单工具
    singles = opt_df[(opt_df['strategy'] == 'single') & (opt_df['virus_type'] == vt_label)]
    all_tools = sorted(singles['subset'].unique())

    for idx, lvl in enumerate(TAX_LEVELS):
        ax = axes[idx]
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_thetagrids(np.degrees(angles), METRIC_LABELS, fontsize=8, weight='bold')
        ax.set_rlabel_position(0)
        ax.set_yticks([20, 40, 60, 80, 100])
        ax.set_yticklabels(["20", "40", "60", "80", "100"], color="grey", fontsize=7)
        ax.set_ylim(0, 100)

        # 所有工具线
        for ti, tool in enumerate(all_tools):
            m = get_metrics_row(opt_df, tool, 'single', lvl, vt_label)
            m += m[:1]
            ax.plot(angles_pad, m, color=TOOL_COLORS[ti % len(TOOL_COLORS)],
                    linewidth=1.5, label=tool[:12] if idx == 0 else "",
                    alpha=0.6)

        # 最优集成方案（加粗虚线）
        if best_df is not None:
            br = best_df[(best_df['level'] == lvl) &
                         (best_df['virus_type'] == vt_label)]
            if not br.empty:
                b = br.iloc[0]
                bm = get_metrics_row(opt_df, b['subset'], b['strategy'], lvl, vt_label)
                bm += bm[:1]
                ax.plot(angles_pad, bm, color='#d62728', linewidth=3.0, linestyle='--',
                        label='Best Ens.' if idx == 0 else "")
                ax.fill(angles_pad, bm, color='#d62728', alpha=0.08)

        ax.set_title(f"{lvl}", fontsize=12, weight='bold', pad=12)
        ax.grid(True, color='#e0e0e0', linestyle='--', linewidth=0.5)

    fig.subplots_adjust(bottom=0.10, hspace=0.30, wspace=0.30)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=min(6, len(handles)),
               fontsize=9, frameon=True)
    plt.suptitle(f"Multi-Metric Radar — {vt_title}", fontsize=16, weight='bold', y=1.01)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [radar] {out_path}")


# ══════════════════════════════════════
# Main
# ══════════════════════════════════════

PALETTE = ['#4C72B0','#55A868','#C44E52','#7E6148','#E64B35','#F39B7F','#8B5CF6','#EC4899','#F59E0B','#06B6D4']


def main():
    parser = argparse.ArgumentParser(description='分类集成结果总结+可视化')
    parser.add_argument('--opt-tsv', required=True, help='ensemble_optimization.tsv')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    opt_df = pd.read_csv(args.opt_tsv, sep='\t')

    # 1. 三指标文本 → 同时生成 best_per_level
    summarize_three_metrics(opt_df, args.outdir)

    # 从刚生成的文件读取 best_df，按 combined_score 排重
    best_path = os.path.join(args.outdir, 'ensemble_best_per_level.tsv')
    if os.path.exists(best_path):
        best_df = pd.read_csv(best_path, sep='\t')
        # 对每个 (level, virus_type) 只保留 sort_metric=='combined_score' 的行，去重
        if 'sort_metric' in best_df.columns:
            best_df = best_df[best_df['sort_metric'] == 'combined_score'].copy()
        best_df = best_df.drop_duplicates(subset=['level', 'virus_type'])
    else:
        best_df = None

    # 2. 条形图
    plot_compare_bars(opt_df, args.outdir)

    # 3. 折线图
    if best_df is not None:
        plot_performance_curves(opt_df, best_df, args.outdir)

    # 4. 假阳性柱状图
    plot_decoy_fpr(opt_df, args.outdir)

    # 5-6. 雷达图
    if best_df is not None:
        draw_radar_grid(opt_df, best_df, 'known', 'Known',
                        os.path.join(args.outdir, 'radar_grid_known_4x2.png'))
        draw_radar_grid(opt_df, best_df, 'novel', 'Novel',
                        os.path.join(args.outdir, 'radar_grid_novel_4x2.png'))

    print(f"\nDone. {args.outdir}")


if __name__ == '__main__':
    main()
