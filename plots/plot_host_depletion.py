#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
===============================================================================
📊 宏病毒组去宿主消融实验可视化与统计工具 (Host Depletion Plotter & Summarizer)
===============================================================================
此脚本用于读取包含宿主与病毒 Reads 统计信息的 TSV 文件：
  1. 自动生成高质量（SCI 级别）的 1x3 组合图 (.pdf 和 .png)
     * 智能适应：自动检测深度梯度数量，单一深度自动切换为分类误差散点图！
  2. 自动聚合并计算去宿主效率，生成核心统计表格 (summary.tsv)
===============================================================================
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--input", type=str, default="host_depletion_detail.tsv",
                        help="输入的 TSV 统计详情文件路径 (默认: host_depletion_detail.tsv)")
    parser.add_argument("-o", "--out_dir", type=str, default="step5_host_free_analysis",
                        help="输出结果的保存目录 (默认: step5_host_free_analysis)")
    parser.add_argument("--prefix", type=str, default="Host_Depletion_Ablation_Study",
                        help="输出图表的文件名前缀 (默认: Host_Depletion_Ablation_Study)")
    parser.add_argument("-p", "--palette", type=str, default="tab10",
                        help="Seaborn 颜色主题 (默认: tab10)")
    parser.add_argument("--dpi", type=int, default=300, help="输出 PNG 的分辨率 (默认: 300)")

    return parser.parse_args()

def main():
    args = parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    if not os.path.exists(args.input):
        print(f"❌ 错误: 找不到输入文件 '{args.input}'。")
        sys.exit(1)

    print(f"⏳ 正在读取数据: {args.input} ...")
    df = pd.read_csv(args.input, sep='\t')
    
    group_order = sorted(df['Group'].dropna().unique().tolist())
    df['Group'] = pd.Categorical(df['Group'], categories=group_order, ordered=True)

    # --- 1. 生成统计汇总表格 ---
    print("⏳ 正在计算去宿主统计汇总表...")
    summary = df.groupby('Group', observed=False)[['Host_Reads', 'Virus_Reads', 'Total_Reads']].sum().reset_index()
    base_row = summary.iloc[0]
    
    summary['Host_Total'] = summary['Host_Reads']
    summary['Virus_Total'] = summary['Virus_Reads']
    summary['Host_Depletion_Rate(%)'] = ((1 - summary['Host_Reads'] / base_row['Host_Reads']) * 100).round(2).clip(lower=0.0)
    summary['Virus_Retention_Rate(%)'] = ((summary['Virus_Reads'] / base_row['Virus_Reads']) * 100).round(2)
    
    base_pct = base_row['Virus_Reads'] / base_row['Total_Reads'] if base_row['Total_Reads'] > 0 else 1.0
    current_pct = summary['Virus_Reads'] / summary['Total_Reads']
    summary['Virus_Enrichment_Fold'] = (current_pct / base_pct).round(1)

    summary_out = os.path.join(args.out_dir, "host_depletion_summary.tsv")
    final_summary = summary[['Group', 'Host_Total', 'Virus_Total', 'Host_Depletion_Rate(%)', 'Virus_Retention_Rate(%)', 'Virus_Enrichment_Fold']]
    final_summary.to_csv(summary_out, sep='\t', index=False)

    # --- 2. 绘制 1x3 组合图 ---
    plt.style.use('default')
    sns.set_theme(style="ticks", context="paper", font_scale=1.2)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    my_palette = sns.color_palette(args.palette, n_colors=len(group_order))

    # A. 宿主清除箱线图
    sns.boxplot(data=df, x='Group', y='Host_Reads', hue='Group', ax=axes[0], palette=my_palette, showfliers=False, legend=False, width=0.5, linecolor='gray', linewidth=1.5)
    sns.stripplot(data=df, x='Group', y='Host_Reads', hue='Group', ax=axes[0], palette=my_palette, size=5, alpha=0.6, jitter=True, legend=False)
    axes[0].set_yscale('log')
    axes[0].set_title('A. Host Reads Depletion (Log Scale)', fontweight='bold')
    axes[0].tick_params(axis='x', rotation=30)
    axes[0].set_xlabel('Depletion Strategy')

    # B. 病毒保留率
    sns.boxplot(data=df, x='Group', y='Virus_Retention', hue='Group', ax=axes[1], palette=my_palette, showfliers=False, legend=False, width=0.5, linecolor='gray', linewidth=1.5)
    sns.stripplot(data=df, x='Group', y='Virus_Retention', hue='Group', ax=axes[1], palette=my_palette, size=5, alpha=0.6, jitter=True, legend=False)
    axes[1].set_ylim(90, 105)
    axes[1].set_title('B. Viral Reads Retention Rate', fontweight='bold')
    axes[1].tick_params(axis='x', rotation=30)
    axes[1].set_xlabel('Depletion Strategy')

    # C. 丰度富集曲线 (🌟 动态自适应核心逻辑)
    actual_depths = sorted(df['LoD_Factor'].unique())
    axes[2].set_yscale('log')
    
    if len(actual_depths) == 1:
        # 只有一个梯度时：采用 分类误差散点图 (Categorical Point Plot)
        single_depth = actual_depths[0]
        
        # 绘制均值与误差棒 (标准差)
        sns.pointplot(data=df, x='Group', y='Virus_Pct', hue='Group', ax=axes[2], 
                      palette=my_palette, errorbar='sd', capsize=0.1, markers="D", 
                      linewidth=2, legend=False)
        # 叠加半透明散点展示底层分布
        sns.stripplot(data=df, x='Group', y='Virus_Pct', hue='Group', ax=axes[2], 
                      palette=my_palette, size=4, alpha=0.4, jitter=True, legend=False)
        
        axes[2].set_title(f'C. Viral Abundance Enrichment (Depth = {single_depth}x)', fontweight='bold')
        axes[2].set_xlabel('Depletion Strategy')
        axes[2].tick_params(axis='x', rotation=30)
        
    else:
        # 有多个梯度时：保持经典的折线图
        sns.lineplot(data=df, x='LoD_Factor', y='Virus_Pct', hue='Group', marker='o', ax=axes[2], palette=my_palette)
        axes[2].set_xscale('log')
        axes[2].set_xticks(actual_depths)
        axes[2].set_xticklabels(actual_depths, rotation=45)
        axes[2].minorticks_off()
        axes[2].set_title('C. Viral Abundance Enrichment', fontweight='bold')
        axes[2].set_xlabel('Target Sequencing Depth (x)')
        axes[2].legend(title='Strategy', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.prefix}.pdf"), dpi=args.dpi, bbox_inches='tight')
    plt.savefig(os.path.join(args.out_dir, f"{args.prefix}.png"), dpi=args.dpi, bbox_inches='tight')
    
    print(f"🎉 绘图完成！图表与汇总表已成功保存至: {args.out_dir}/")

if __name__ == "__main__":
    main()
