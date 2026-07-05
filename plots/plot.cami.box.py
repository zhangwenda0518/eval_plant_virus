#!/usr/bin/env python3
"""
【宏基因组评估绘图工具 - 顶刊 2x2 布局版 (含全局图例)】
提取 OPAL results.tsv 核心指标，绘制 Nature Methods 风格箱线图。
支持在无 GUI 界面的 Linux 服务器上静默运行，并支持布局自定义。
"""

import os
import argparse
import sys
import pandas as pd

# 强制 Matplotlib 使用非交互式后台（解决服务器 xcb 崩溃问题）
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches # 用于手动构建完美的全局图例

# ==========================================
# 核心配置：指标映射字典
# ==========================================
METRIC_MAPPING = {
    "Completeness": "Sensitivity",
    "Purity": "Precision",
    "False positives": "FalsePositives",
    "L1 norm error": "L1 Error"
}
METRIC_ORDER = ["Sensitivity", "Precision", "FalsePositives", "L1 Error"]

def plot_panel_b(input_tsv, output_fig, target_rank, layout="2x2", custom_tools=None):
    if not os.path.exists(input_tsv):
        print(f"[错误] 找不到输入文件: {input_tsv}")
        sys.exit(1)

    print(f"📥 正在读取数据: {input_tsv}")
    df = pd.read_csv(input_tsv, sep='\t')

    # 过滤数据
    df = df[(df['rank'] == target_rank) & (df['tool'] != 'Gold standard')]

    if df.empty:
        print(f"[错误] 在数据中未找到 rank 为 '{target_rank}' 的预测结果。")
        print("这通常是因为你的原始数据中没有这一层级。请检查文件或尝试 --rank species")
        sys.exit(1)

    # 重命名指标
    df = df[df['metric'].isin(METRIC_MAPPING.keys())].copy()
    df['metric'] = df['metric'].map(METRIC_MAPPING)

    # 确定软件排序
    available_tools = df['tool'].unique().tolist()
    if custom_tools:
        tools_to_plot = [t for t in custom_tools if t in available_tools]
    else:
        tools_to_plot = sorted(available_tools)

    if not tools_to_plot:
        print("[错误] 没有有效的软件用于绘图。")
        sys.exit(1)

    df = df[df['tool'].isin(tools_to_plot)]

    # 开始绘图
    print(f"🎨 正在后台静默绘制箱线图 (级别: {target_rank}, 布局: {layout})...")
    sns.set_theme(style="ticks", font_scale=1.1)

    palette = sns.color_palette("Set2", len(tools_to_plot))

    # 根据用户选择的布局创建画布
    if layout == "2x2":
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8), sharex=True)
        axes_flat = axes.flatten() 
    else: # 1x4
        fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(8, 12), sharex=True)
        axes_flat = axes 

    for i, metric in enumerate(METRIC_ORDER):
        ax = axes_flat[i]
        metric_data = df[df['metric'] == metric]

        # 画箱线图 (保持 legend=False，我们稍后在外部统一加)
        sns.boxplot(data=metric_data, x='tool', y='value', order=tools_to_plot,
                    hue='tool', hue_order=tools_to_plot, palette=palette,
                    ax=ax, width=0.5, fliersize=0, boxprops=dict(alpha=0.6), legend=False)

        # 画散点抖动
        sns.stripplot(data=metric_data, x='tool', y='value', order=tools_to_plot,
                      hue='tool', hue_order=tools_to_plot, palette=palette,
                      ax=ax, size=5, jitter=0.2, dodge=False, linewidth=1,
                      edgecolor='auto', legend=False)

        ax.set_title(metric, fontweight='bold', fontsize=14)
        ax.set_ylabel("Value", fontsize=12)
        ax.set_xlabel("")
        sns.despine(ax=ax)

        # 修复 Warning：先设定 xticks，再设定 xticklabels
        if layout == "2x2":
            if i >= 2: 
                ax.set_xticks(range(len(tools_to_plot)))
                ax.set_xticklabels(tools_to_plot, rotation=45, ha='right', fontweight='bold')
            else: 
                ax.tick_params(labelbottom=False)
        else: # 1x4 布局
            if i == 3: 
                ax.set_xticks(range(len(tools_to_plot)))
                ax.set_xticklabels(tools_to_plot, rotation=45, ha='right', fontweight='bold')
            else:
                ax.tick_params(labelbottom=False)

    # 🌟 新增：手动构建全局图例 🌟
    # 根据已有的 palette 和 tools 创建图例色块
    legend_handles = [mpatches.Patch(color=palette[j], label=tool) for j, tool in enumerate(tools_to_plot)]
    
    # 将图例放置在整个画布的中心偏右侧外围
    fig.legend(handles=legend_handles, title="Tool", title_fontsize=13, fontsize=12,
               loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    # 调整子图间距并添加大标题
    plt.tight_layout(rect=[0, 0, 1, 0.96]) 
    fig.suptitle(f'Profiler Performance at {target_rank.capitalize()} Level', fontsize=18, fontweight='bold')

    # 保存图片 (bbox_inches='tight' 确保外围的图例不会被裁剪掉)
    fig.savefig(output_fig, dpi=300, bbox_inches='tight')
    print(f"✅ 图表已成功保存至: {os.path.abspath(output_fig)}")

def main():
    parser = argparse.ArgumentParser(
        description="提取 OPAL 结果并绘制顶刊风格箱线图 (默认 2x2 布局)。",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input', required=True, help="OPAL 的 results.tsv 文件路径")
    parser.add_argument('-o', '--output', default="Performance_Boxplots.png", help="输出文件名")
    parser.add_argument('-r', '--rank', default="species", choices=["phylum", "class", "order", "family", "genus", "species", "strain"])
    parser.add_argument('-t', '--tools', nargs='+', help="指定软件顺序")
    parser.add_argument('-l', '--layout', default="2x2", choices=["2x2", "1x4"], help="图片布局方式")

    args = parser.parse_args()
    plot_panel_b(args.input, args.output, args.rank, args.layout, args.tools)

if __name__ == "__main__":
    main()
