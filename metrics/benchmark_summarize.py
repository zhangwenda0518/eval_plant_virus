#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Virus Assembly Benchmark — Phase 4: 汇总 & 可视化 (定制进阶版)
================================================================================

功能:
  读取前三个 Phase 的 TSV 输出 → 合并 → 绘图 + 排名表 + 加权得分表
  🌟 智能读取 selected_viruses.tsv，自动进行 Segmented / Non-Segmented 分组对比
  🌟 深度折线图启用极具学术感的 Log-Scale 坐标轴

用法:
  python scripts/benchmark_summarize.py -d benchmark_results --mode 7
"""

import os, sys, re, argparse, warnings

import pandas as pd
import numpy as np
warnings.filterwarnings('ignore', category=RuntimeWarning)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# ── 非对称误差棒：下限不低于0 ──
def _safe_errorbar(ax, x, y, std, label, color, marker='o', **kw):
    lower = np.minimum(y, std)
    upper = std
    ax.errorbar(x, y, yerr=[lower, upper], marker=marker, capsize=3.5,
                label=label, color=color, linewidth=1.5, markersize=5, alpha=0.85, **kw)

# ── 统一图例样式：双行三列，增大间距防重叠 ──
LEGEND_DEFAULTS = dict(
    loc='upper center',
    ncol=3,
    frameon=True,
    facecolor='white',
    edgecolor='#e2e2e2',
    framealpha=0.9,
    fontsize=9.5,
    columnspacing=1.8,
    handletextpad=0.6,
)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import (
    MY_PALETTE, METAQUAST_METRICS, SCORE_WEIGHTS, BIG_BETTER,
    MUT_ORDER, DEPTH_ORDER, FILE_PATTERNS
)

PLOT_METRICS = [
    "Genome fraction (%)",
    "Duplication ratio",
    "# mismatches per 100 kbp",
    "# indels per 100 kbp",
    "# misassemblies",
    "Largest_alignment_relative",
    "NGA50_relative",
]

# 🌟 坐标系刻度格式化函数 (1000000 -> 1M)
def depth_formatter(x, pos):
    if x >= 1000000: return f'{int(x/1000000)}M'
    if x >= 1000: return f'{int(x/1000)}K'
    return str(int(x))

def load_data(data_dir):
    """加载三个 Phase 的 TSV 文件"""
    mq_file = os.path.join(data_dir, "metaquast_summary.tsv")
    chim_file = os.path.join(data_dir, "chimeric_summary.tsv")
    res_file = os.path.join(data_dir, "resource_summary.tsv")

    df_mq = None; df_chim = None; df_res = None
    missing = []

    if os.path.exists(mq_file):
        df_mq = pd.read_csv(mq_file, sep='\t')
        print(f"✅ metaquast_summary.tsv: {len(df_mq)} 行")
    else:
        missing.append("metaquast_summary.tsv")

    if os.path.exists(chim_file):
        df_chim = pd.read_csv(chim_file, sep='\t')
        print(f"✅ chimeric_summary.tsv: {len(df_chim)} 行")
    else:
        missing.append("chimeric_summary.tsv")

    if os.path.exists(res_file):
        df_res = pd.read_csv(res_file, sep='\t')
        print(f"✅ resource_summary.tsv: {len(df_res)} 行")
    else:
        missing.append("resource_summary.tsv")

    if missing:
        print(f"⚠️  缺少: {', '.join(missing)}，相关图表将跳过")
    return df_mq, df_chim, df_res


def prepare_data(df_mq):
    """准备病毒级分析数据：长度归一化 → 3重复平均 → 返回 (per_virus_mean_df, per_virus_metrics)
    每行 = 一个病毒在一个深度下一个工具的 3rep 均值"""
    df = df_mq.copy()
    # 转数值
    num_cols = list(PLOT_METRICS) + ["Reference length", "Largest alignment", "NGA50"]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

    # 计算长度归一化指标（病毒级，聚合前）
    if "Reference length" in df.columns:
        rl = df["Reference length"].fillna(1)
        if "Largest alignment" in df.columns:
            df["Largest_alignment_relative"] = (df["Largest alignment"].fillna(0) / rl * 100).clip(0, 100)
        if "NGA50" in df.columns:
            df["NGA50_relative"] = (df["NGA50"].fillna(0) / rl * 100).clip(0, 100)

    # 按 (Tool, Depth[, Background_Ratio], Virus) 聚合 3 重复 → mean, std
    gcols = ['Tool', 'Depth', 'Virus']
    if 'Background_Ratio' in df.columns and df['Background_Ratio'].nunique(dropna=True) > 1:
        gcols.insert(2, 'Background_Ratio')
    mcols = [c for c in PLOT_METRICS if c in df.columns]
    # 保留 Reference length 用于指标表（不参与绘图）
    if "Reference length" in df.columns:
        mcols.append("Reference length")
    if not mcols:
        return pd.DataFrame(), []
    agg = df.groupby(gcols)[mcols].agg(['mean', 'std']).reset_index()
    agg.columns = [f"{a}_{b}" if b else a for a, b in agg.columns]
    # 恢复 mean 列为原名
    rename = {}
    for m in mcols:
        rename[f"{m}_mean"] = m
        rename[f"{m}_std"] = f"{m}_sd"
    agg.rename(columns=rename, inplace=True)

    # 只保留 mean 列 + Tool/Depth/Virus，作为主数据；sd 列用于后续
    data_cols = gcols + mcols
    sd_cols = [f"{m}_sd" for m in mcols]
    per_virus = agg[data_cols].copy()
    per_virus_sd = agg[gcols + sd_cols].copy()
    print(f"✅ 病毒级数据: {len(per_virus)} 行 ({per_virus['Virus'].nunique()} viruses)")
    return per_virus, mcols, per_virus_sd


def merge_summary(df_mq, df_chim, df_res):
    """合并三个数据源为统一的 benchmark_summary，并自动打上节段标签"""
    merged = df_mq.copy() if df_mq is not None else pd.DataFrame()

    # 🌟 修改点 1：读取 reference TSV，自动识别 Segmented 病毒
    if not merged.empty and 'Virus' in merged.columns:
        seg_map = {}
        ref_tsv = "step1_eval_viruses/selected_viruses.tsv"
        if os.path.exists(ref_tsv):
            try:
                ref_df = pd.read_csv(ref_tsv, sep='\t')
                for _, row in ref_df.iterrows():
                    acc = str(row['accession']).split('.')[0]
                    gtype = str(row['genome_type'])
                    if "Segmented" in gtype or "Segmented_RNA" in gtype:
                        seg_map[acc] = "Segmented"
                    else:
                        seg_map[acc] = "Non-Segmented"
            except Exception as e:
                print(f"⚠️ 读取 {ref_tsv} 失败 (但这不影响继续运行): {e}")
        
        def assign_group(v):
            v_clean = str(v).split('.')[0]
            return seg_map.get(v_clean, "Non-Segmented")
        
        merged['Group'] = merged['Virus'].apply(assign_group)

    # 🌟 修改点 2：合并嵌合数据 (去掉 Group 作为键，因为嵌合表是样本维度的)
    if df_chim is not None and not df_chim.empty:
        chim_cols = ['Tool', 'Mutation_Rate', 'Depth', 'Rep',
                     'Total_Contigs', 'Chimeric_Count', 'Chimeric_Rate_pct']
        chim_sel = df_chim[[c for c in chim_cols if c in df_chim.columns]]
        if not merged.empty:
            merged = merged.merge(chim_sel, on=['Tool', 'Mutation_Rate', 'Depth', 'Rep'], how='left')
        else:
            merged = chim_sel

    # 🌟 修改点 3：合并资源数据 (同理去掉 Group)
    if df_res is not None and not df_res.empty:
        res_cols = [c for c in ['Tool', 'Mutation_Rate', 'Depth', 'Rep',
                                'Wall_Time_s', 'Max_RSS_MB', 'CPU_pct',
                                'IO_In_MB', 'IO_Out_MB'] if c in df_res.columns]
        res_sel = df_res[res_cols]
        if not merged.empty:
            merged = merged.merge(res_sel, on=['Tool', 'Mutation_Rate', 'Depth', 'Rep'], how='left')
        else:
            merged = res_sel

    return merged


def _get_facet_var(df, candidates, varname):
    unique_vals = df[varname].dropna().unique()
    available = [v for v in candidates if v in set(unique_vals)]
    if not available:
        available = sorted(unique_vals)
    return available if len(available) > 1 else None


def _plot_box_facet(df, col, x_var, row_var, row_vals, active_tools, out_dir, title_suffix):
    """通用箱线图 (分面)"""
    data = df.dropna(subset=[col])
    data = data[data['Tool'].isin(active_tools)]
    if data.empty: return

    if row_vals and len(row_vals) > 1:
        g = sns.catplot(
            data=data, x=x_var, y=col, hue="Tool", row=row_var, kind="box",
            row_order=row_vals, hue_order=active_tools,
            palette=MY_PALETTE, height=3.2, aspect=1.6, fliersize=1.2, legend_out=True
        )
        g.fig.suptitle(f"{col}  (row = {row_var})", fontsize=16, y=1.02, fontweight='bold')
    else:
        g = sns.catplot(
            data=data, x=x_var, y=col, hue="Tool", kind="box",
            hue_order=active_tools, palette=MY_PALETTE,
            height=4, aspect=1.8, fliersize=1.2, legend_out=True
        )
        g.fig.suptitle(f"{col}  ({x_var} only)", fontsize=16, y=1.02, fontweight='bold')
    
    g.set_axis_labels(x_var, col)
    # 对于 Depth 轴，稍微转义下让其更好看
    if x_var == "Depth":
        for ax in g.axes.flat:
            labels = ax.get_xticklabels()
            new_labels = []
            for l in labels:
                try:
                    v = float(l.get_text())
                    new_labels.append(depth_formatter(v, 0))
                except:
                    new_labels.append(l.get_text())
            ax.set_xticklabels(new_labels)

    safe = col.replace("/", "per").replace(" ", "_").replace("#", "Num").replace("(", "").replace(")", "").replace("%", "pct")
    fname = os.path.join(out_dir, f"Plot_{x_var}_{row_var}_{safe}.png")
    g.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   📊 {os.path.basename(fname)}")


def plot_quality_metrics(df, metrics, active_tools, out_dir):
    mut_vals = _get_facet_var(df, MUT_ORDER, "Mutation_Rate")
    depth_vals = _get_facet_var(df, DEPTH_ORDER, "Depth")

    for col in metrics:
        if col not in df.columns: continue
        if depth_vals and depth_vals != [0]:
            if mut_vals:
                _plot_box_facet(df, col, "Depth", "Mutation_Rate", mut_vals, active_tools, out_dir, "row=MutationRate")
            else:
                _plot_box_facet(df, col, "Depth", None, None, active_tools, out_dir, "depth-only")

        if mut_vals and depth_vals and depth_vals != [0]:
            df2 = df.copy()
            df2['Mutation_Rate'] = pd.Categorical(df2['Mutation_Rate'], categories=MUT_ORDER, ordered=True)
            _plot_box_facet(df2, col, "Mutation_Rate", "Depth", depth_vals, active_tools, out_dir, "row=Depth")


def plot_chimeric_both(df, active_tools, out_dir):
    """嵌合率折线图（X轴微调错位避免重叠）"""
    if 'Chimeric_Rate_pct' not in df.columns: return
    xvar = 'Background_Ratio' if ('Background_Ratio' in df.columns and df['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    data = df[df['Tool'].isin(active_tools)].dropna(subset=['Chimeric_Rate_pct', xvar])
    if data.empty: return
    fig, ax = plt.subplots(figsize=(9, 5.5))
    present_tools = [t for t in active_tools if t in data['Tool'].unique()]
    n = len(present_tools)
    for i, tool in enumerate(present_tools):
        td = data[data['Tool'] == tool]
        if td.empty: continue
        agg = td.groupby(xvar)['Chimeric_Rate_pct'].agg(['mean', 'std']).reset_index()
        valid = agg[xvar].notna()
        if xvar == 'Depth': valid = valid & (agg[xvar] > 0)
        agg = agg[valid].sort_values(xvar)
        if agg.empty: continue
        offset = 1.05 ** (i - (n - 1) / 2)
        x_shifted = agg[xvar] * offset
        lower = np.minimum(agg['mean'], agg['std'])
        ax.errorbar(x_shifted, agg['mean'], yerr=[lower, agg['std']], marker='o',
                    capsize=3.5, label=tool, color=MY_PALETTE.get(tool, '#888'),
                    linewidth=1.5, markersize=5, alpha=0.85)
    ax.set_ylim(bottom=0)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.set_xlabel(xvar.replace("_", " ").title(), fontsize=11, labelpad=10)
    ax.set_ylabel("Chimeric Rate (%)", fontsize=11, labelpad=10)
    ax.set_title("Chimeric Contig Rate", fontsize=12, fontweight='bold', pad=12)
    ax.grid(True, alpha=0.15, linestyle='--')
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        fig.legend(handles[:len(present_tools)], labels[:len(present_tools)],
                  bbox_to_anchor=(0.5, -0.05), **LEGEND_DEFAULTS)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(os.path.join(out_dir, "Plot_Chimeric_Line.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   📊 Plot_Chimeric_Line.png")

    # ── 方案 B：分组条形图 ──
    fig2, ax2 = plt.subplots(figsize=(9, 5.5))
    data[xvar] = data[xvar].astype(str)
    x_order_str = [str(x) for x in sorted(df[xvar].dropna().unique(), key=lambda x: float(x))]
    sns.barplot(data=data, x=xvar, y='Chimeric_Rate_pct', hue='Tool',
                order=x_order_str, hue_order=present_tools, palette=MY_PALETTE,
                edgecolor='black', linewidth=0.7, capsize=0.1,
                err_kws={'linewidth': 1.1}, ax=ax2)
    if ax2.legend_ is not None:
        ax2.legend_.remove()
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel(xvar.replace("_", " ").title(), fontsize=11, labelpad=10)
    ax2.set_ylabel("Chimeric Rate (%)", fontsize=11, labelpad=10)
    ax2.set_title("Chimeric Contig Rate (Grouped Bar)", fontsize=12, fontweight='bold', pad=12)
    ax2.grid(True, alpha=0.15, axis='y', linestyle='--')
    handles2, labels2 = ax2.get_legend_handles_labels()
    if handles2:
        fig2.legend(handles2[:len(present_tools)], labels2[:len(present_tools)],
                   bbox_to_anchor=(0.5, -0.05), **LEGEND_DEFAULTS)
    fig2.tight_layout(rect=[0, 0.05, 1, 1])
    fig2.savefig(os.path.join(out_dir, "Plot_Chimeric_Bar.png"), dpi=300, bbox_inches='tight')
    plt.close(fig2)
    print(f"   📊 Plot_Chimeric_Bar.png")


def plot_resource(df_res, active_tools, out_dir):
    """资源消耗图：1×3 一行三列（时间/内存/CPU）"""
    if df_res is None or df_res.empty: return
    xvar = 'Background_Ratio' if ('Background_Ratio' in df_res.columns and df_res['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    specs = [('Wall_Time_s', 'Wall Time (s)'), ('Max_RSS_MB', 'Max RSS (MB)'), ('CPU_pct', 'CPU (%)')]
    specs = [(m, yl) for m, yl in specs if m in df_res.columns]
    if not specs: return

    fig, axes = plt.subplots(1, len(specs), figsize=(5.5 * len(specs), 5))
    axes = [axes] if len(specs) == 1 else axes
    handles, labels = [], []

    for ax, (metric, ylabel) in zip(axes, specs):
        data = df_res[df_res['Tool'].isin(active_tools)].dropna(subset=[metric, xvar])
        if data.empty: continue
        agg = data.groupby(['Tool', xvar])[metric].agg(['mean', 'std']).reset_index()
        for tool in active_tools:
            sub = agg[agg['Tool'] == tool].sort_values(xvar)
            if sub.empty: continue
            line = ax.errorbar(sub[xvar], sub['mean'], yerr=sub['std'],
                       marker='o', capsize=4, label=tool,
                       color=MY_PALETTE.get(tool, '#888'),
                       linewidth=1.5, markersize=5)
            if not handles:
                handles, labels = ax.get_legend_handles_labels()
        if xvar == 'Depth':
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(depth_formatter))
        ax.set_xlabel(xvar.replace("_", " ").title(), fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.grid(True, alpha=0.3)

    if handles:
        fig.legend(handles, labels, loc='lower center', fontsize=8, ncol=len(active_tools),
                  frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.suptitle("Resource Consumption", fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    fig.savefig(os.path.join(out_dir, "Plot_Resource.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   📊 Plot_Resource.png")


def plot_radar(df, active_tools, out_dir, df_res=None):
    """雷达图：含组装质量 7 指标 + 运行时间 + 内存（越小越好，取倒数）"""
    radar_metrics = [c for c in PLOT_METRICS if c in df.columns]
    if not radar_metrics: return
    agg = df.groupby('Tool')[radar_metrics].mean().reset_index()
    # 成功率（GF≥50% 病毒占比）
    if 'Virus' in df.columns and 'Genome fraction (%)' in df.columns:
        total_v = df['Virus'].nunique()
        succ = df.groupby('Tool').apply(
            lambda g: g[g['Genome fraction (%)'] >= 50]['Virus'].nunique() / total_v * 100
        ).reset_index(name='Success Rate (%)')
        agg = agg.merge(succ, on='Tool', how='left')
        radar_metrics.append('Success Rate (%)')
    # 合并资源指标
    if df_res is not None and not df_res.empty:
        res_agg = df_res.groupby('Tool')[['Wall_Time_s', 'Max_RSS_MB']].mean().reset_index()
        agg = agg.merge(res_agg, on='Tool', how='left')
        for mc in ['Wall_Time_s', 'Max_RSS_MB']:
            if mc in agg.columns:
                radar_metrics.append(mc)
    agg = agg[agg['Tool'].isin(active_tools)]
    scaled = agg.copy()
    smaller_better = {"# mismatches per 100 kbp", "# indels per 100 kbp", "# misassemblies",
                      "Duplication ratio", "Wall_Time_s", "Max_RSS_MB"}
    for col in radar_metrics:
        mn, mx = scaled[col].min(), scaled[col].max()
        if mx == mn: scaled[col] = 0.5
        elif col in BIG_BETTER or col == 'Success Rate (%)': scaled[col] = (scaled[col] - mn) / (mx - mn)
        elif col in smaller_better: scaled[col] = 1.0 - (scaled[col] - mn) / (mx - mn)
        else: scaled[col] = 0.5

    tools_list = scaled['Tool'].tolist()
    n = len(radar_metrics)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(projection='polar'))
    for i, tool in enumerate(tools_list):
        vals = scaled.iloc[i][radar_metrics].tolist()
        vals += vals[:1]
        color = MY_PALETTE.get(tool, '#888')
        ax.fill(angles, vals, alpha=0.01, color=color)
        ax.plot(angles, vals, 'o-', linewidth=2, label=tool, color=color, markersize=5)

    label_map = {"Wall_Time_s": "1 / Runtime", "Max_RSS_MB": "1 / Memory",
                 "Largest_alignment_relative": "Largest Align / RefLen",
                 "NGA50_relative": "NGA50 / RefLen",
                 "Success Rate (%)": "Success Rate"}
    labels = []
    for m in radar_metrics:
        if m in label_map: labels.append(label_map[m])
        else: labels.append(m.replace(" (%)","").replace("# ","").replace(" per 100 kbp","/100kbp").replace("_"," ").title())
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.grid(True, color='#e0e0e0', linestyle='--')
    ax.set_title("Radar: Assembly Quality & Resource Efficiency\n(Normalized: Outward is Better)",
                 fontsize=13, fontweight='bold', pad=30)
    ax.legend(loc='upper left', bbox_to_anchor=(1.1, 1.05), fontsize=9, frameon=True, edgecolor='#e2e2e2')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "Plot_Radar.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   📊 Plot_Radar.png")


def plot_heatmap(df, active_tools, out_dir):
    """2×2 四合一热图：GF | SuccessRate / N50Rel | LgRel"""
    if 'Genome fraction (%)' not in df.columns: return
    data = df.dropna(subset=['Genome fraction (%)'])
    data = data[data['Tool'].isin(active_tools)]
    colvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    gb_cols = ['Tool', colvar]

    def _pivot(data, gb, col, metric, fmt='.1f', vmax=100):
        p = data.groupby(gb)[metric].mean().reset_index()
        mat = p.pivot(index='Tool', columns=col, values=metric)
        return mat.reindex(index=[t for t in active_tools if t in mat.index])

    def _pivot_success(data, gb, col, total_v):
        p = data.groupby(gb).apply(
            lambda g: g[g['Genome fraction (%)'] >= 50]['Virus'].nunique() / total_v * 100
        ).reset_index(name='Success_Rate')
        mat = p.pivot(index='Tool', columns=col, values='Success_Rate')
        return mat.reindex(index=[t for t in active_tools if t in mat.index])

    total_v = data['Virus'].nunique()
    specs = [
        ('Genome fraction (%)', 'Genome Fraction (%)', '.1f', 100),
        ('Success_Rate',        'Success Rate (GF≥50%)', '.0f', 100),
        ('NGA50_relative',      'NGA50 / RefLen (%)',    '.1f', 100),
        ('Largest_alignment_relative', 'Largest Align / RefLen (%)', '.1f', 100),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    for ax, (metric, title, fmt, vmax) in zip(axes.flat, specs):
        if metric == 'Success_Rate':
            mat = _pivot_success(data, gb_cols, colvar, total_v)
        elif metric in data.columns:
            mat = _pivot(data, gb_cols, colvar, metric)
        else:
            ax.set_visible(False)
            continue
        if mat.empty:
            ax.set_visible(False)
            continue
        sns.heatmap(mat, annot=True, fmt=fmt, cmap='YlOrRd', ax=ax,
                    vmin=0, vmax=vmax, linewidths=0.5, cbar_kws={'shrink': 0.8})
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel(colvar.replace("_", " ").title(), fontsize=9)
        ax.set_ylabel("")

    fig.suptitle(f"Heatmap: Assembly Quality — Tool × {colvar}", fontsize=15, fontweight='bold', y=1.01)
    fig.tight_layout()
    fpath = os.path.join(out_dir, "Plot_Heatmap_2x2.png")
    fig.savefig(fpath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   📊 Plot_Heatmap_2x2.png")


# ── 指标显示名和标签 ──
_METRIC_LABELS = {
    "Genome fraction (%)":        ("Genome Fraction (%)",       "GF"),
    "Duplication ratio":          ("Duplication Ratio",         "Dup"),
    "# mismatches per 100 kbp":   ("Mismatches / 100kbp",       "Mism"),
    "# indels per 100 kbp":       ("Indels / 100kbp",           "Indel"),
    "# misassemblies":            ("Misassemblies",              "MisAsm"),
    "Largest_alignment_relative": ("Largest Align / RefLen (%)", "LgRel"),
    "NGA50_relative":             ("NGA50 / RefLen (%)",         "N50Rel"),
}


def plot_metric_box(data, active_tools, out_dir):
    """箱线图：legend=False禁内部图例，统一底端图例，计数指标Y轴下限=0"""
    if 'Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1:
        xvar = 'Background_Ratio'
        x_order = sorted(data[xvar].dropna().unique(), key=lambda x: float(x))
    else:
        xvar = 'Depth'
        x_order = sorted([d for d in data[xvar].unique() if d > 0])

    for metric, (ylabel, tag) in _METRIC_LABELS.items():
        if metric not in data.columns: continue
        sub = data[data['Tool'].isin(active_tools)].dropna(subset=[metric]).copy()
        if sub.empty: continue
        sub[xvar] = sub[xvar].astype(str)
        present_x = [str(x) for x in x_order if str(x) in sub[xvar].unique()]
        present_tools = [t for t in active_tools if t in sub['Tool'].unique()]
        if not present_tools or not present_x: continue

        fig, ax = plt.subplots(figsize=(9, 6))
        sns.boxplot(data=sub, x=xvar, y=metric, hue='Tool', order=present_x,
                    hue_order=present_tools, palette=MY_PALETTE,
                    showfliers=False, linewidth=1.1, ax=ax, legend=False)
        sns.stripplot(data=sub, x=xvar, y=metric, hue='Tool', order=present_x,
                      hue_order=present_tools, palette=MY_PALETTE,
                      dodge=True, alpha=0.2, size=2.5, jitter=0.15, ax=ax, legend=False)
        if any(kw in ylabel for kw in ["Mismatches", "Indels", "Misassemblies", "Ratio", "Chimeric"]):
            ax.set_ylim(bottom=0)
        ax.set_xlabel(xvar.replace("_", " ").title(), fontsize=11, labelpad=10)
        ax.set_ylabel(ylabel, fontsize=11, labelpad=10)
        ax.set_title(f"{ylabel} (n={sub['Virus'].nunique()} viruses, 3-rep mean)", fontsize=12, fontweight='bold', pad=12)
        ax.grid(True, alpha=0.15, axis='y', linestyle='--')
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            fig.legend(handles[:len(present_tools)], labels[:len(present_tools)],
                      bbox_to_anchor=(0.5, -0.04), **LEGEND_DEFAULTS)
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        fig.savefig(os.path.join(out_dir, f"box_{tag}.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)
    print(f"   📊 metric box: {len(_METRIC_LABELS)} plots")


def plot_sd_table(data, sd_df, active_tools, out_dir):
    """输出三重复 SD 汇总表"""
    xvar = 'Background_Ratio' if ('Background_Ratio' in sd_df.columns and sd_df['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    rows = []
    for tool in active_tools:
        td = sd_df[sd_df['Tool'] == tool]
        for xv in sorted(td[xvar].dropna().unique()):
            dd = td[td[xvar] == xv]
            row = {'Tool': tool, xvar: xv}
            for metric in _METRIC_LABELS:
                sc = f"{metric}_sd"
                if sc in dd.columns:
                    row[f"{_METRIC_LABELS[metric][1]}_meanSD"] = dd[sc].mean()
                    row[f"{_METRIC_LABELS[metric][1]}_medianSD"] = dd[sc].median()
            rows.append(row)
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(out_dir, "sd_summary.tsv"), sep='\t', index=False)
        print(f"   📊 sd_summary.tsv")


def plot_metric_matrix(data, active_tools, out_dir):
    """3行×N列矩阵面板：row=GF/LgRel/N50Rel, col=Tool, X=Depth(log)或Background_Ratio"""
    metrics_spec = [
        ("Genome fraction (%)",        "GF (%)"),
        ("Largest_alignment_relative",  "Largest Align / RefLen (%)"),
        ("NGA50_relative",              "NGA50 / RefLen (%)"),
    ]
    xvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    present_tools = [t for t in active_tools if t in data['Tool'].unique()]
    if not present_tools: return

    num_rows, num_cols = 3, len(present_tools)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(2.4 * num_cols, 8.5),
                             sharex='col', sharey='row')
    if num_rows == 1: axes = np.atleast_2d(axes)
    elif num_cols == 1: axes = np.atleast_2d(axes).T

    for row_idx, (metric_col, ylabel) in enumerate(metrics_spec):
        for col_idx, tool in enumerate(present_tools):
            ax = axes[row_idx, col_idx]
            td = data[(data['Tool'] == tool) & data[metric_col].notna()]
            if td.empty: continue
            agg = td.groupby(xvar)[metric_col].agg(['mean', 'std']).reset_index()
            if xvar == 'Depth': agg = agg[agg[xvar] > 0]
            agg = agg.dropna(subset=[xvar]).sort_values(xvar)
            if agg.empty: continue
            _safe_errorbar(ax, agg[xvar], agg['mean'], agg['std'], tool, MY_PALETTE.get(tool, '#888'))
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(depth_formatter) if xvar == 'Depth' else ticker.FormatStrFormatter('%g'))
            ax.set_ylim(0, 105)
            ax.grid(True, alpha=0.15, linestyle='--')
            ax.tick_params(axis='both', which='major', labelsize=9)
            if row_idx == 0:
                ax.set_title(tool, fontsize=11, fontweight='bold', pad=12)
            if col_idx == 0:
                ax.set_ylabel(ylabel, fontsize=10.5, fontweight='bold', labelpad=8)
            if row_idx == num_rows - 1:
                ax.set_xlabel(xvar.replace("_", " ").title(), fontsize=9.5, labelpad=8)
                plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    fig.suptitle("Core Metrics Matrix (GF, LgRel, N50Rel) — Tool × Metric", fontsize=14, fontweight='bold', y=0.97)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(out_dir, "Panel_Metrics_vs_Tools.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   📊 Panel_Metrics_vs_Tools.png")


def plot_tool_panel(data, active_tools, out_dir):
    """一张大图：1×7 子图（7指标），每个子图含所有工具线"""
    all_metrics = list(_METRIC_LABELS.items())
    xvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    fig, axes = plt.subplots(2, 4, figsize=(18, 9.5))
    axes_flat = axes.flatten()

    for idx, (metric, (ylabel, tag)) in enumerate(all_metrics):
        ax = axes_flat[idx]
        if metric not in data.columns:
            ax.set_visible(False); continue
        sub = data.dropna(subset=[metric])
        for tool in active_tools:
            td = sub[sub['Tool'] == tool]
            if td.empty: continue
            agg = td.groupby(xvar)[metric].agg(['mean', 'std']).reset_index()
            valid = agg[xvar].notna()
            if xvar == 'Depth': valid = valid & (agg[xvar] > 0)
            agg = agg[valid].sort_values(xvar)
            if agg.empty: continue
            _safe_errorbar(ax, agg[xvar], agg['mean'], agg['std'], tool, MY_PALETTE.get(tool, '#888'))
        ax.set_ylim(bottom=0)
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
        ax.set_title(ylabel, fontsize=10.5, fontweight='bold', pad=8)
        ax.set_xlabel(xvar.replace("_", " ").title(), fontsize=9)
        ax.grid(True, alpha=0.15, linestyle='--')
        ax.tick_params(axis='both', which='major', labelsize=8.5)

    if len(all_metrics) < 8:
        axes_flat[-1].set_visible(False)
    # 图例提取移到所有工具绘制之后（修复只显示1个工具的Bug）
    handles, labels = axes_flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.02), **LEGEND_DEFAULTS)
    fig.suptitle("Assembly Quality Panel (n=96 viruses, 3-rep mean)", fontsize=14, fontweight='bold', y=0.97)
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    fig.savefig(os.path.join(out_dir, "Panel_All_Tools.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   📊 Panel_All_Tools.png")


def sanitize(s):
    return s.replace("/", "per").replace(" ", "_").replace("#", "Num").replace("(", "").replace(")", "").replace("%", "pct")

_RADAR_METRICS = list(_METRIC_LABELS.keys())


def _norm_for_radar(series, bigger_better):
    """归一化到 [0,1]：越大越好则 min-max，越小越好则 1/(x+1) 倒数后 min-max"""
    v = series.fillna(0)
    if bigger_better:
        mn, mx = v.min(), v.max()
        return pd.Series(0.5, index=v.index) if mx == mn else (v - mn) / (mx - mn)
    inv = 1.0 / (v.clip(lower=0) + 1)
    mn, mx = inv.min(), inv.max()
    return pd.Series(0.5, index=inv.index) if mx == mn else (inv - mn) / (mx - mn)


SEGMENT_MARKERS = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'X', 'P', 'd', '8']


def _build_seg_map():
    """从 selected_viruses.tsv 构建 Segmented_RNA 的 species→[accessions] 映射"""
    ref_tsv = "step1_eval_viruses/selected_viruses.tsv"
    seg2sp = {}
    sp2segs = {}
    if not os.path.exists(ref_tsv):
        return seg2sp, sp2segs
    try:
        ref = pd.read_csv(ref_tsv, sep='\t')
        for _, row in ref.iterrows():
            if str(row.get('genome_type', '')) != 'Segmented_RNA':
                continue
            acc = str(row['accession']).split('.')[0]
            sp = str(row['species'])
            seg2sp[acc] = sp
            sp2segs.setdefault(sp, []).append(acc)
    except Exception:
        pass
    return seg2sp, sp2segs


def plot_per_virus_panel(data, active_tools, out_dir, viruses_filter=None):
    """--viruses 指定：单病毒3×3面板（7指标×Depth折线）+ 雷达图。节段病毒自动展开同species全部segment"""
    if 'Virus' not in data.columns:
        return
    all_metrics = list(_METRIC_LABELS.items())
    ncols = 3; nrows = 3
    xvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'

    all_viruses_set = set(data['Virus'].unique())
    seg2sp, sp2segs = _build_seg_map()

    # 匹配用户指定的病毒（自动去版本号：NC_002030.1 → NC_002030）
    all_viruses = sorted(all_viruses_set)
    if viruses_filter:
        filters_expanded = set(viruses_filter)
        for vf in viruses_filter:
            base = vf.rsplit('.', 1)[0] if '.' in vf and vf.rsplit('.', 1)[1].isdigit() else None
            if base and base != vf:
                filters_expanded.add(base)
        matched = [v for v in all_viruses if any(
            v.startswith(f) or f in v for f in filters_expanded)]
        # 展开节段病毒 + 去重
        species_groups = {}
        single_viruses = set()
        processed_sp = set()
        for v in matched:
            sp = seg2sp.get(v)
            if sp and sp not in processed_sp:
                siblings = [s for s in sp2segs.get(sp, []) if s in all_viruses_set]
                if len(siblings) > 1:
                    species_groups[sp] = siblings
                    processed_sp.add(sp)
                else:
                    single_viruses.add(v)
            elif not sp:
                single_viruses.add(v)
        all_viruses = sorted(single_viruses)
    else:
        species_groups = {}
        all_viruses = sorted(all_viruses_set)

    if not all_viruses and not species_groups:
        return

    # ── 非节段病毒 ──
    for virus in all_viruses:
        vd = data[data['Virus'] == virus]
        if vd.empty:
            continue

        # ── 3×3 面板（7指标）──
        fig, axes = plt.subplots(nrows, ncols, figsize=(16, 11))
        axes = np.atleast_2d(axes)
        for idx, (metric, (ylabel, tag)) in enumerate(all_metrics):
            r, c = idx // ncols, idx % ncols
            ax = axes[r, c]
            if metric not in vd.columns:
                ax.set_visible(False); continue
            for tool in active_tools:
                td = vd[(vd['Tool'] == tool)].dropna(subset=[metric])
                if td.empty: continue
                agg = td.groupby(xvar)[metric].agg(['mean', 'std']).reset_index()
                agg = (agg[agg[xvar] > 0] if xvar == 'Depth' else agg).sort_values(xvar)
                if agg.empty: continue
                ax.errorbar(agg[xvar], agg['mean'], yerr=agg['std'],
                           marker='o', capsize=3, label=tool,
                           color=MY_PALETTE.get(tool, '#888'),
                           linewidth=1.5, markersize=4, alpha=0.85)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(depth_formatter))
            ax.set_title(ylabel, fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
        for idx in range(len(all_metrics), nrows * ncols):
            axes[idx // ncols, idx % ncols].set_visible(False)
        for c in range(ncols):
            axes[-1, c].set_xlabel(xvar.replace("_", " ").title(), fontsize=9)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc='lower center', fontsize=7,
                      ncol=len(active_tools), frameon=False, bbox_to_anchor=(0.5, -0.02))
        fig.suptitle(virus, fontsize=14, fontweight='bold')
        fig.tight_layout(rect=[0, 0.04, 1, 0.96])
        sv = virus.replace("/", "_").replace(" ", "_")
        fig.savefig(os.path.join(out_dir, f"panel_{sv}.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)

        # ── 雷达图 ──
        available = [(m, lab, bb) for m, (lab, _), bb in zip(
            _RADAR_METRICS,
            [(l, t) for l, t in _METRIC_LABELS.values()],
            [True, False, False, False, False, True, True]) if m in vd.columns]
        vm = vd.groupby('Tool')[list({a[0] for a in available})].mean()
        tools = [t for t in active_tools if t in vm.index]
        if len(tools) < 2:
            continue
        normed = pd.DataFrame({'Tool': tools})
        for metric, label, bigger in available:
            normed[label] = _norm_for_radar(vm[metric].reindex(tools), bigger).values
        n = len(available)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist() + [0]
        fig2, ax2 = plt.subplots(figsize=(7, 7), subplot_kw=dict(projection='polar'))
        for _, row in normed.iterrows():
            vals = [row[l] for _, l, _ in available] + [row[available[0][1]]]
            color = MY_PALETTE.get(row['Tool'], '#888')
            ax2.fill(angles, vals, alpha=0.06, color=color)
            ax2.plot(angles, vals, 'o-', linewidth=2, label=row['Tool'], color=color, markersize=4)
        ax2.set_xticks(angles[:-1])
        ax2.set_xticklabels([lab for _, lab, _ in available], fontsize=8)
        ax2.set_ylim(0, 1.05)
        ax2.set_title(virus, fontsize=12, fontweight='bold', pad=20)
        ax2.legend(loc='upper right', bbox_to_anchor=(1.3, 1.05), fontsize=7)
        fig2.tight_layout()
        fig2.savefig(os.path.join(out_dir, f"radar_{sv}.png"), dpi=200, bbox_inches='tight')
        plt.close(fig2)

    # ── 节段病毒（同 species 多 segment 合并为一张面板）──
    for sp, segments in species_groups.items():
        segs_sorted = sorted(segments)
        vd = data[data['Virus'].isin(segs_sorted)]
        if vd.empty: continue

        # 3×3 面板（多 segment：工具区分颜色，节段区分标记）
        fig, axes = plt.subplots(nrows, ncols, figsize=(16, 11))
        axes = np.atleast_2d(axes)
        for idx, (metric, (ylabel, tag)) in enumerate(all_metrics):
            r, c = idx // ncols, idx % ncols
            ax = axes[r, c]
            if metric not in vd.columns:
                ax.set_visible(False); continue
            for tool in active_tools:
                for si, seg in enumerate(segs_sorted):
                    td = vd[(vd['Tool'] == tool) & (vd['Virus'] == seg)].dropna(subset=[metric])
                    if td.empty: continue
                    agg = td.groupby(xvar)[metric].agg(['mean', 'std']).reset_index()
                    agg = (agg[agg[xvar] > 0] if xvar == 'Depth' else agg).sort_values(xvar)
                    if agg.empty: continue
                    marker = SEGMENT_MARKERS[si % len(SEGMENT_MARKERS)]
                    ax.errorbar(agg[xvar], agg['mean'], yerr=agg['std'],
                               marker=marker, capsize=3,
                               color=MY_PALETTE.get(tool, '#888'),
                               linewidth=1.3, markersize=4, alpha=0.8)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(depth_formatter))
            ax.set_title(ylabel, fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
        for idx in range(len(all_metrics), nrows * ncols):
            axes[idx // ncols, idx % ncols].set_visible(False)
        for c in range(ncols):
            axes[-1, c].set_xlabel(xvar.replace("_", " ").title(), fontsize=9)

        # 双图例：工具=颜色，节段=标记
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        tool_patches = [Patch(color=MY_PALETTE.get(t, '#888'), label=t)
                        for t in active_tools if t in vd['Tool'].unique()]
        seg_handles = [Line2D([0],[0], marker=SEGMENT_MARKERS[i % len(SEGMENT_MARKERS)],
                       color='k', linestyle='None', markersize=6, label=s)
                       for i, s in enumerate(segs_sorted)]
        leg1 = fig.legend(handles=tool_patches, title="Tools",
                         bbox_to_anchor=(0.02, 0), loc='lower left',
                         fontsize=7, title_fontsize=8, ncol=min(4, len(tool_patches)), frameon=False)
        fig.add_artist(leg1)
        fig.legend(handles=seg_handles, title="Segments",
                  bbox_to_anchor=(0.98, 0), loc='lower right',
                  fontsize=7, title_fontsize=8, ncol=min(3, len(segs_sorted)), frameon=False)
        fig.suptitle(f"Segmented | {sp} ({len(segs_sorted)} segments)", fontsize=14, fontweight='bold')
        fig.tight_layout(rect=[0, 0.06, 1, 0.96])
        sv = sp.replace("/", "_").replace(" ", "_")
        fig.savefig(os.path.join(out_dir, f"panel_{sv}.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)

        # 雷达图（取各segment均值）
        available = [(m, lab, bb) for m, (lab, _), bb in zip(
            _RADAR_METRICS,
            [(l, t) for l, t in _METRIC_LABELS.values()],
            [True, False, False, False, False, True, True]) if m in vd.columns]
        vm = vd.groupby('Tool')[list({a[0] for a in available})].mean() if available else pd.DataFrame()
        tools = [t for t in active_tools if t in vm.index] if not vm.empty else []
        if len(tools) >= 2:
            normed = pd.DataFrame({'Tool': tools})
            for metric, label, bigger in available:
                normed[label] = _norm_for_radar(vm[metric].reindex(tools), bigger).values
            n = len(available)
            angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist() + [0]
            fig3, ax3 = plt.subplots(figsize=(7, 7), subplot_kw=dict(projection='polar'))
            for _, row in normed.iterrows():
                vals = [row[l] for _, l, _ in available] + [row[available[0][1]]]
                color = MY_PALETTE.get(row['Tool'], '#888')
                ax3.fill(angles, vals, alpha=0.06, color=color)
                ax3.plot(angles, vals, 'o-', linewidth=2, label=row['Tool'], color=color, markersize=4)
            ax3.set_xticks(angles[:-1])
            ax3.set_xticklabels([lab for _, lab, _ in available], fontsize=8)
            ax3.set_ylim(0, 1.05)
            ax3.set_title(f"{sp} ({len(segs_sorted)} segs)", fontsize=12, fontweight='bold', pad=20)
            ax3.legend(loc='upper right', bbox_to_anchor=(1.3, 1.05), fontsize=7)
            fig3.tight_layout()
            fig3.savefig(os.path.join(out_dir, f"radar_{sv}.png"), dpi=200, bbox_inches='tight')
            plt.close(fig3)

    print(f"   📊 per_virus: {len(all_viruses)} non-seg + {len(species_groups)} segmented")


def plot_virus_metric_table(data, active_tools, out_dir, viruses_filter=None):
    """--viruses 指定时，输出 MetaQUAST 参考对比指标表（分病毒分深度，论文 Table 格式）"""
    if not viruses_filter or data is None or data.empty:
        return

    filters_expanded = set(viruses_filter)
    for vf in viruses_filter:
        base = vf.rsplit('.', 1)[0] if '.' in vf and vf.rsplit('.', 1)[1].isdigit() else None
        if base:
            filters_expanded.add(base)
    matched = [v for v in sorted(data['Virus'].unique()) if any(
        v.startswith(f) or f in v for f in filters_expanded)]
    if not matched:
        print("   ⚠️ 未匹配到指定病毒"); return

    # 指标顺序
    metric_spec = [
        ("Genome fraction (%)",        "Genome fraction",    "%.3f"),
        ("Largest_alignment_relative",  "LgAlign/RefLen(%)",  "%.1f"),
        ("NGA50_relative",              "NGA50/RefLen(%)",    "%.1f"),
        ("Duplication ratio",           "Dup ratio",          "%.2f"),
        ("# mismatches per 100 kbp",    "Mism/100kbp",        "%.1f"),
        ("# indels per 100 kbp",        "Indel/100kbp",       "%.1f"),
        ("# misassemblies",             "MisAsm",             "%.1f"),
    ]
    metric_cols = [m for m, _, _ in metric_spec if m in data.columns]

    for virus in matched:
        vd = data[data['Virus'] == virus]
        xvar = 'Background_Ratio' if ('Background_Ratio' in vd.columns and vd['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
        xvals = sorted(vd[xvar].dropna().unique())
        if not xvals: continue

        ref_len = int(vd['Reference length'].mean()) if 'Reference length' in vd.columns else None

        lines = []
        xlabel = xvar.replace("_", " ").title()
        header = xlabel + "\t" + "\t".join(f"{lab}" for _, lab, _ in metric_spec if _[0] in metric_cols)
        lines.append(f"# Virus: {virus}" + (f"  (RefLen={ref_len}bp)" if ref_len else ""))
        lines.append(header)

        for xv in xvals:
            dd = vd[vd[xvar] == xv]
            depth_printed = False
            for tool in active_tools:
                td = dd[dd['Tool'] == tool]
                if td.empty:
                    continue
                vals = []
                for m, _, fmt in metric_spec:
                    if m not in metric_cols:
                        continue
                    v = td[m].mean()
                    vals.append(fmt % v if not np.isnan(v) else "-")
                tag = f"{xv}×" if xvar == 'Depth' and not depth_printed else (f"bg{xv}" if xvar == 'Background_Ratio' and not depth_printed else "")
                lines.append(f"{tag}\t{tool}\t" + "\t".join(vals))
                depth_printed = True

        out_path = os.path.join(out_dir, f"table_{sanitize(virus)}.tsv")
        with open(out_path, 'w') as f:
            f.write("\n".join(lines) + "\n")
        print(f"   📊 {os.path.basename(out_path)}")
    print(f"   📊 virus_metrics_table: {len(matched)} viruses")


def plot_assembly_success_rate(data, active_tools, out_dir, threshold=80.0):
    """成功率折线图：Y=GF≥threshold%的病毒占比，消除0%和100%双极化"""
    if 'Genome fraction (%)' not in data.columns: return
    xvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    sub = data[data['Tool'].isin(active_tools)]
    records = []
    for (tool, xv), grp in sub.groupby(['Tool', xvar]):
        n_virus = grp['Virus'].nunique()
        if n_virus == 0: continue
        n_success = grp[grp['Genome fraction (%)'] >= threshold]['Virus'].nunique()
        records.append({'Tool': tool, xvar: xv, 'Success_Rate': n_success / n_virus * 100})
    df = pd.DataFrame(records)
    if df.empty: return
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for tool in [t for t in active_tools if t in df['Tool'].unique()]:
        td = df[df['Tool'] == tool].sort_values(xvar)
        ax.plot(td[xvar], td['Success_Rate'], marker='o', linewidth=2, markersize=6,
                label=tool, color=MY_PALETTE.get(tool, '#888'))
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(depth_formatter) if xvar == 'Depth' else ticker.FormatStrFormatter('%g'))
    ax.set_ylim(-5, 105)
    ax.set_xlabel("Sequencing Depth" if xvar == 'Depth' else "Background Ratio", fontsize=11, labelpad=8)
    ax.set_ylabel(f"Genomes Assembled (GF ≥ {int(threshold)}%)", fontsize=11, labelpad=8)
    ax.set_title(f"Assembly Success Rate (n={sub['Virus'].nunique()} viruses)", fontsize=12, fontweight='bold', pad=12)
    ax.grid(True, alpha=0.15, linestyle='--')
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, bbox_to_anchor=(0.5, -0.05), **LEGEND_DEFAULTS)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(os.path.join(out_dir, f"Plot_Success_Rate_GF{int(threshold)}.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   📊 Plot_Success_Rate_GF{int(threshold)}.png")


def plot_metric_violin(data, active_tools, out_dir):
    """小提琴图代替箱线图：展示0%和100%双峰密度分布"""
    xvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    for metric, (ylabel, tag) in _METRIC_LABELS.items():
        if metric not in data.columns: continue
        sub = data[data['Tool'].isin(active_tools)].dropna(subset=[metric]).copy()
        if sub.empty: continue
        sub[xvar] = sub[xvar].astype(str)
        x_order_str = [str(x) for x in sorted(sub[xvar].dropna().unique(), key=lambda x: float(x))]
        present_tools = [t for t in active_tools if t in sub['Tool'].unique()]
        fig, ax = plt.subplots(figsize=(9.5, 6))
        sns.violinplot(data=sub, x=xvar, y=metric, hue='Tool', order=x_order_str,
                       hue_order=present_tools, palette=MY_PALETTE,
                       inner='quartile', linewidth=1.0, ax=ax, bw_adjust=0.4, legend=False)
        if "Ratio" in ylabel or "fraction" in ylabel.lower() or "Align" in ylabel:
            ax.set_ylim(-5, 105)
        elif any(kw in ylabel for kw in ["Mismatches", "Indels", "Misassemblies"]):
            ax.set_ylim(bottom=-0.5)
        ax.set_xlabel("Depth" if xvar == 'Depth' else "Background Ratio", fontsize=11, labelpad=8)
        ax.set_ylabel(ylabel, fontsize=11, labelpad=8)
        ax.set_title(f"{ylabel} Density (n={sub['Virus'].nunique()} viruses)", fontsize=12, fontweight='bold', pad=12)
        ax.grid(True, alpha=0.15, axis='y', linestyle='--')
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            fig.legend(handles[:len(present_tools)], labels[:len(present_tools)],
                      bbox_to_anchor=(0.5, -0.05), **LEGEND_DEFAULTS)
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        fig.savefig(os.path.join(out_dir, f"violin_{tag}.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)
    print(f"   📊 violin: {len(_METRIC_LABELS)} plots")


def plot_virus_count_by_gf_thresholds(data, active_tools, out_dir):
    """分组条形图：X=Depth(或BgRatio), 分组=Tool, Y=达标病毒数, 多张(GF30/50/70/90/100)"""
    if 'Genome fraction (%)' not in data.columns or 'Virus' not in data.columns: return
    xvar = 'Background_Ratio' if ('Background_Ratio' in data.columns and data['Background_Ratio'].nunique(dropna=True) > 1) else 'Depth'
    total = data['Virus'].nunique()
    present_tools = [t for t in active_tools if t in data['Tool'].unique()]
    x_labels = [depth_formatter(x, 0) if xvar == 'Depth' else str(x)
                for x in sorted(data[xvar].dropna().unique())]

    for threshold in [30, 50, 70, 90, 100]:
        records = []
        for tool in present_tools:
            td = data[data['Tool'] == tool]
            for xv in sorted(td[xvar].dropna().unique()):
                sub = td[td[xvar] == xv]
                n = sub[sub['Genome fraction (%)'] >= threshold]['Virus'].nunique()
                label = depth_formatter(xv, 0) if xvar == 'Depth' else str(xv)
                records.append({'Tool': tool, 'XV': label, 'Count': n})
        dfp = pd.DataFrame(records)
        if dfp.empty: continue

        fig, ax = plt.subplots(figsize=(max(9, len(x_labels) * 1.6), 6))
        sns.barplot(data=dfp, x='XV', y='Count', hue='Tool', order=x_labels,
                    hue_order=present_tools, palette=MY_PALETTE,
                    edgecolor='black', linewidth=0.6, ax=ax)
        # 移除 Seaborn 自动生成的右侧内部图例
        if ax.legend_ is not None:
            ax.legend_.remove()
        # 柱顶标注数字
        for p in ax.patches:
            h = p.get_height()
            if h > 0:
                ax.annotate(f'{int(h)}', (p.get_x() + p.get_width() / 2., h),
                           ha='center', va='bottom', fontsize=7, color='#333',
                           fontweight='bold', xytext=(0, 2), textcoords='offset points')
        ax.set_ylim(0, total * 1.15)
        ax.set_xlabel(xvar.replace("_", " ").title(), fontsize=11, labelpad=8)
        ax.set_ylabel(f"Viruses (GF ≥ {threshold}%)", fontsize=11, labelpad=8)
        ax.set_title(f"Assembly Count at GF≥{threshold}% (total: {total} viruses)", fontsize=12, fontweight='bold', pad=12)
        ax.grid(True, alpha=0.1, axis='y', linestyle='--')
        # 统一底端左对齐一行图例
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            fig.legend(handles[:len(present_tools)], labels[:len(present_tools)],
                      loc='lower left', bbox_to_anchor=(0.02, -0.02),
                      ncol=len(present_tools), frameon=True, facecolor='white',
                      edgecolor='#e2e2e2', fontsize=8, columnspacing=1.2, handletextpad=0.5)
        fig.tight_layout(rect=[0, 0.05, 1, 1])
        fig.savefig(os.path.join(out_dir, f"count_GF{threshold}_vs_{xvar}.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

        # 输出病毒名单（组装成功 vs 未组装）
        list_rows = []
        for xv in sorted(data[xvar].dropna().unique()):
            for tool in present_tools:
                sub = data[(data['Tool'] == tool) & (data[xvar] == xv)]
                assembled = sub[sub['Genome fraction (%)'] >= threshold]['Virus'].unique()
                unassembled = set(sub['Virus'].unique()) - set(assembled)
                list_rows.append({
                    xvar: xv, 'Tool': tool, f'GF≥{threshold}%_Assembled': ','.join(sorted(assembled)),
                    f'GF≥{threshold}%_Unassembled': ','.join(sorted(unassembled)),
                })
        pd.DataFrame(list_rows).to_csv(os.path.join(out_dir, f"virus_list_GF{threshold}_{xvar}.tsv"),
                                       sep='\t', index=False)
    print(f"   📊 count_GFthreshold: 5 plots + 5 virus lists")


def compute_ranking_score(df, active_tools, out_dir):
    """成功率门控得分：GF≥50%病毒占比 × 质量分（低覆盖工具无高分）"""
    if 'Genome fraction (%)' not in df.columns or 'Virus' not in df.columns: return
    gf_key = 'Genome fraction (%)'
    quality_keys = [c for c in SCORE_WEIGHTS if c in df.columns and c != gf_key]
    if not quality_keys: return

    # 成功率：每个工具在所有条件下 GF≥50% 的病毒占比
    total_viruses = df['Virus'].nunique()
    success_df = df.groupby('Tool').apply(
        lambda g: g[g[gf_key] >= 50]['Virus'].nunique() / total_viruses
    ).reset_index(name='Success_Rate')

    # 质量分
    num_cols = [c for c in quality_keys if c in df.columns]
    df_num = df[num_cols + ['Tool']].copy()
    for c in num_cols:
        df_num[c] = pd.to_numeric(df_num[c], errors='coerce')
    qual_agg = df_num.groupby('Tool')[num_cols].mean().reset_index()

    merged = success_df.merge(qual_agg, on='Tool')
    merged = merged[merged['Tool'].isin(active_tools)]
    if merged.empty: return

    scaled = merged.copy()
    for col in quality_keys:
        mn, mx = scaled[col].min(), scaled[col].max()
        if mx == mn: scaled[col] = 0.5
        elif col in BIG_BETTER: scaled[col] = (scaled[col] - mn) / (mx - mn)
        else: scaled[col] = (mx - scaled[col]) / (mx - mn)

    raw_w = [SCORE_WEIGHTS.get(m, 0.1) for m in quality_keys]
    total_w = sum(raw_w)
    quality_w = [w / total_w for w in raw_w]
    scaled['Quality_Score'] = scaled[quality_keys].dot(quality_w)
    scaled['Weighted_Score'] = scaled['Success_Rate'] * scaled['Quality_Score'] * 10

    score_table = scaled[['Tool', 'Success_Rate', 'Quality_Score', 'Weighted_Score']].sort_values(
        'Weighted_Score', ascending=False)
    score_table.to_csv(os.path.join(out_dir, "benchmark_weighted_score.tsv"), sep='\t', index=False, float_format='%.4f')
    print(f"   📊 benchmark_weighted_score.tsv")
    print(f"\n  🏆 加权综合得分 (Success% × Quality, 0-10):")
    print(f"    {'Tool':20s}  {'Success%':>8s}  {'Quality':>7s}  {'Score':>7s}")
    for _, r in score_table.iterrows():
        print(f"    {r['Tool']:20s}  {r['Success_Rate']:.4f}  {r['Quality_Score']:.4f}  {r['Weighted_Score']:.2f}")


def main():
    parser = argparse.ArgumentParser(description="🧬 Phase 4: 汇总 & 可视化")
    parser.add_argument("-d", "--data-dir", default="benchmark_results")
    parser.add_argument("--mode", choices=['6', '10'], default='6',
                        help="6=基础6工具 10=全部含merge")
    parser.add_argument("--phase", choices=['all', 'assembly', 'chimeric', 'resource', 'ranking'],
                        default='all', help="只跑指定模块")
    parser.add_argument("--viruses", nargs="*", default=None,
                        help="指定病毒 accession（空格分隔），生成单病毒面板+雷达图+指标表")
    args = parser.parse_args()

    do_all  = args.phase in ('all', 'assembly')
    do_chim = args.phase in ('all', 'chimeric')
    do_res  = args.phase in ('all', 'resource')
    do_rank = args.phase in ('all', 'ranking')

    os.makedirs(args.data_dir, exist_ok=True)
    plots_dir = os.path.join(args.data_dir, "plots")
    box_dir   = os.path.join(plots_dir, "01_metric_box")
    panel_dir = os.path.join(plots_dir, "02_panel")
    virus_dir = os.path.join(plots_dir, "03_per_virus")
    for d in [box_dir, panel_dir, virus_dir]:
        os.makedirs(d, exist_ok=True)

    if args.mode == '10':
        active_tools = list(FILE_PATTERNS.keys())
    else:  # mode '6' (default)
        active_tools = ["Megahit", "RNAViralSPAdes", "Penguin", "MetaViralSPAdes", "Trinity", "RNABloom"]

    print("📥 加载数据...")
    df_mq, df_chim, df_res = load_data(args.data_dir)

    if do_all:
        print("🔗 准备病毒级数据...")
        virus_data, mcols, sd_data = prepare_data(df_mq)
        print("🔗 合并样本级数据（嵌合/资源）...")
        summary = merge_summary(df_mq, df_chim, df_res)
        for col in PLOT_METRICS + list(_METRIC_LABELS.keys()) + ['Wall_Time_s', 'Max_RSS_MB', 'CPU_pct', 'Chimeric_Rate_pct']:
            if col in summary.columns:
                summary[col] = pd.to_numeric(summary[col], errors='coerce')
    else:
        virus_data, mcols, sd_data = None, [], None
        if do_chim or do_res:
            summary = merge_summary(df_mq, df_chim, df_res)
            for col in PLOT_METRICS + ['Wall_Time_s', 'Max_RSS_MB', 'CPU_pct', 'Chimeric_Rate_pct']:
                if col in summary.columns:
                    summary[col] = pd.to_numeric(summary[col], errors='coerce')

    print(f"\n🎨 生成图表...")
    sns.set_theme(style="ticks", font_scale=1.05)

    if do_all:
        print("\n  === 1. Violin ===")
        plot_metric_violin(virus_data, active_tools, box_dir)
        print("\n  === 2. Success Rate ===")
        plot_assembly_success_rate(virus_data, active_tools, box_dir, threshold=80.0)
        plot_assembly_success_rate(virus_data, active_tools, box_dir, threshold=50.0)
        print("\n  === 3. SD Table ===")
        plot_sd_table(virus_data, sd_data, active_tools, args.data_dir)
        print("\n  === 4. Tool Panel ===")
        plot_tool_panel(virus_data, active_tools, panel_dir)
        print("\n  === 5. Matrix Panel ===")
        plot_metric_matrix(virus_data, active_tools, panel_dir)
        print("\n  === GF Threshold Curves ===")
        plot_virus_count_by_gf_thresholds(virus_data, active_tools, box_dir)
        if args.viruses:
            print("\n  === Per-Virus ===")
            plot_per_virus_panel(virus_data, active_tools, virus_dir, args.viruses)
            plot_virus_metric_table(virus_data, active_tools, virus_dir, args.viruses)
        print("\n  === Global Radar & Heatmap ===")
        plot_radar(virus_data, active_tools, plots_dir, df_res)
        plot_heatmap(virus_data, active_tools, plots_dir)

    if do_chim and df_chim is not None and not df_chim.empty:
        print("\n  === Chimeric Rate ===")
        plot_chimeric_both(summary, active_tools, plots_dir)

    if do_res and df_res is not None and not df_res.empty:
        print("\n  === Resource ===")
        plot_resource(df_res, active_tools, plots_dir)

    if do_rank and virus_data is not None:
        print("\n  === Ranking ===")
        compute_ranking_score(virus_data, active_tools, args.data_dir)

    print(f"\n🎉 完成！图表目录: {plots_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
