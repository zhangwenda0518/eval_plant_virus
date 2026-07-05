#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="⚡ 宏病毒组计算资源消耗可视化与汇总工具 (Boxplot + Stripplot)")
    parser.add_argument("-i", "--log_dir", type=str, default="step5_host_free_logs", help="日志根目录")
    parser.add_argument("-o", "--out_dir", type=str, default="step5_host_free_analysis", help="输出目录")
    parser.add_argument("--prefix", type=str, default="Resource_Usage_Benchmarking", help="输出前缀")
    parser.add_argument("-p", "--palette", type=str, default="tab10", help="颜色主题")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. 搜索并读取日志
    search_pattern = os.path.join(args.log_dir, "*/host_depletion_resource_usage.tsv")
    log_files = glob.glob(search_pattern)
    
    if not log_files:
        print(f"❌ 错误: 未找到资源日志文件。")
        sys.exit(1)

    all_data = []
    for f in log_files:
        group_name = os.path.basename(os.path.dirname(f))
        try:
            df_tmp = pd.read_csv(f, sep='\t')
            df_tmp['Group'] = group_name
            all_data.append(df_tmp)
        except Exception as e:
            pass

    df = pd.concat(all_data, ignore_index=True)
    group_order = sorted(df['Group'].unique().tolist())
    df['Group'] = pd.Categorical(df['Group'], categories=group_order, ordered=True)

    # --- 2. 生成汇总表 (新增样本数量统计) ---
    res_summary = df.groupby('Group', observed=False).agg(
        Sample_Count=('Sample', 'count'),          # 统计读入了多少个样本
        Time_Mean_Sec=('Elapsed_Seconds', 'mean'), # 平均耗时
        Time_Std=('Elapsed_Seconds', 'std'),       # 耗时波动 (标准差)
        Mem_Peak_MB=('Peak_Mem_MB', 'max'),        # 取内存的最高峰值
        CPU_Avg=('Avg_CPUs', 'mean')               # 平均并行度
    ).round(2).reset_index()
    res_summary.to_csv(os.path.join(args.out_dir, "host_depletion_resource_summary.tsv"), sep='\t', index=False)

    # --- 3. 绘制 1x2 分布对比图 (箱线图 + 散点图) ---
    plt.style.use('default')
    sns.set_theme(style="ticks", context="paper", font_scale=1.2)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    my_palette = sns.color_palette(args.palette, n_colors=len(group_order))

    # A. 运行时间 (恢复为 Boxplot)
    sns.boxplot(data=df, x='Group', y='Elapsed_Seconds', hue='Group', ax=axes[0], 
                palette=my_palette, showfliers=False, legend=False, 
                width=0.5, linecolor='gray', linewidth=1.5)
    # 添加半透明散点，展示每一个样本的具体位置
    sns.stripplot(data=df, x='Group', y='Elapsed_Seconds', hue='Group', ax=axes[0], 
                  palette=my_palette, size=4, alpha=0.6, jitter=True, legend=False)
    axes[0].set_title('A. Elapsed Time per Sample', fontweight='bold')
    axes[0].set_ylabel('Time (Seconds)')
    axes[0].tick_params(axis='x', rotation=30)

    # B. 内存峰值 (恢复为 Boxplot)
    sns.boxplot(data=df, x='Group', y='Peak_Mem_MB', hue='Group', ax=axes[1], 
                palette=my_palette, showfliers=False, legend=False, 
                width=0.5, linecolor='gray', linewidth=1.5)
    sns.stripplot(data=df, x='Group', y='Peak_Mem_MB', hue='Group', ax=axes[1], 
                  palette=my_palette, size=4, alpha=0.6, jitter=True, legend=False)
    axes[1].set_title('B. Peak Memory Usage', fontweight='bold')
    axes[1].set_ylabel('Memory (MB)')
    axes[1].tick_params(axis='x', rotation=30)

    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.prefix}.pdf"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(args.out_dir, f"{args.prefix}.png"), dpi=300, bbox_inches='tight')
    
    print(f"🎉 资源消耗分布图绘制完成！图表与汇总表已保存至: {args.out_dir}/")

if __name__ == "__main__":
    main()
