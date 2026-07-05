#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
【宏基因组评估神器】多软件资源消耗 (时间与内存) 全自动绘图工具
=============================================================================
特性 (Features):
  1. 全自动递归搜索目录下的资源消耗表 (默认: sample_resource_usage.tsv)
  2. 自动清洗并提取模拟数据集条件 (如: abundance, mut0~mut30)
  3. 智能算法分类映射 (Mapper, Aligner, Kmer)，并在图中以不同几何形状区分
  4. 一键生成科研级双轴并排对比图，内存轴自动采用 Log10 对数缩放
  5. 底部居中生成美观的全局独立图例

版本 (Version): 1.1.0 (终极完美发行版)
=============================================================================
"""

__version__ = "1.1.0"

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="📊 【宏基因组评估】多软件资源消耗自动绘图工具 (算法分类+形状映射版)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
使用示例:
  1. 搜索当前目录 (默认):
     plot_identify.resources.py
  
  2. 指定输入目录和输出前缀:
     plot_identify.resources.py -i step7_identify -o step7_identify_analysis/my_eval
"""
    )
    parser.add_argument('-i', '--input', default='.', help="输入根目录 (默认: 当前目录 '.')")
    parser.add_argument('-f', '--filename', default='sample_resource_usage.tsv', help="目标文件名")
    parser.add_argument('-o', '--out_prefix', default='resource_evaluation', help="输出前缀 (可包含路径)")
    parser.add_argument('--palette', default='tab20', help="线条颜色方案 (默认: tab20)")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args()

def main():
    args = parse_args()
    search_dir = Path(args.input)

    # 如果输出路径包含文件夹，自动创建该文件夹
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"===========================================================")
    print(f"🚀 宏基因组比对工具资源评估分析 (v{__version__})")
    print(f"===========================================================")
    print(f"🔍 正在目录 '{search_dir.resolve()}' 中搜索 '{args.filename}' ...")
    
    found_files = list(search_dir.rglob(args.filename))

    if not found_files:
        print(f"❌ 搜索完毕：未找到任何目标文件。")
        sys.exit(1)

    dfs = []
    print(f"✅ 共找到 {len(found_files)} 个数据文件:")
    for file_path in found_files:
        try:
            dfs.append(pd.read_csv(file_path, sep='\t'))
            print(f"  [成功] {file_path.relative_to(search_dir) if search_dir.is_absolute() else file_path}")
        except Exception as e:
            print(f"  [错误] 无法解析 {file_path}: {e}")

    if not dfs:
        sys.exit(1)

    df = pd.concat(dfs, ignore_index=True)

    # ==========================================
    # 🌟 定义算法大类映射字典
    # ==========================================
    algo_mapping = {
        'KALLISTO': 'Mapper', 'SALMON': 'Mapper',
        'BOWTIE2': 'Aligner', 'BWA': 'Aligner', 'MINIMAP2': 'Aligner',
        'STROBEALIGN': 'Aligner', 'BWA-MEM2': 'Aligner', 'HISAT2': 'Aligner',
        'KRAKEN2': 'Kmer', 'KRAKENUNIQ': 'Kmer', 'CENTRIFUGER': 'Kmer',
        'KUNPENG': 'Kmer', 'KAIJU': 'Kmer', 'KRAKEN2X': 'Kmer',
        'METABULI': 'Kmer', 'GANON': 'Kmer', 'SYLPH': 'Kmer'
    }
    
    df['Category'] = df['Tool'].str.upper().map(algo_mapping).fillna('Other')
    marker_dict = {'Mapper': 'o', 'Aligner': 's', 'Kmer': 'D', 'Other': 'X'}

    # ==========================================
    # 🌟 数据清洗与排序
    # ==========================================
    df['Condition'] = df['Sample'].str.replace('eval2b_', '').str.replace('eval2a_', '').str.replace('_PE', '')
    order = ['abundance', 'mut0', 'mut5', 'mut10', 'mut15', 'mut30']
    extra_conditions = [c for c in df['Condition'].unique() if c not in order]
    final_order = order + sorted(extra_conditions)
    
    df['Condition'] = pd.Categorical(df['Condition'], categories=final_order, ordered=True)
    df = df.sort_values(['Condition', 'Tool'])

    unique_tools = df['Tool'].unique()
    unique_conditions = df['Condition'].dropna().unique()

    # --- 终端报告基本信息 ---
    print(f"\n📊 数据合并完毕！包含以下基本信息:")
    print(f"  - 识别到 {len(unique_tools)} 种工具: {', '.join(unique_tools)}")
    print(f"  - 评估条件横轴: {', '.join(unique_conditions)}")
    
    # 找出极限消耗任务
    max_time_idx = df['Time(s)'].idxmax()
    max_mem_idx = df['Peak_Memory(MB)'].idxmax()
    print(f"  - 最耗时记录: {df.loc[max_time_idx, 'Tool']} @ {df.loc[max_time_idx, 'Condition']} ({df.loc[max_time_idx, 'Time(s)']} 秒)")
    print(f"  - 最高内存记录: {df.loc[max_mem_idx, 'Tool']} @ {df.loc[max_mem_idx, 'Condition']} ({df.loc[max_mem_idx, 'Peak_Memory(MB)']} MB)")
    print(f"-----------------------------------------------------------")

    # ==========================================
    # 🌟 开始绘图
    # ==========================================
    print("🎨 正在生成并渲染评估曲线图...")
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # ---- 图 A：运行时间 ----
    sns.lineplot(
        data=df, x='Condition', y='Time(s)', 
        hue='Tool', style='Category', markers=marker_dict, dashes=False, 
        markersize=9, linewidth=2.5, palette=args.palette, ax=axes[0]
    )
    axes[0].set_title('Execution Time across Different Datasets', fontsize=15, fontweight='bold')
    axes[0].set_xlabel('Dataset (Mutation Rate / Abundance)', fontsize=13)
    axes[0].set_ylabel('Time (Seconds)', fontsize=13)
    axes[0].tick_params(axis='x', rotation=30)
    if axes[0].get_legend(): axes[0].get_legend().remove()

    # ---- 图 B：峰值内存 (Log 坐标) ----
    sns.lineplot(
        data=df, x='Condition', y='Peak_Memory(MB)', 
        hue='Tool', style='Category', markers=marker_dict, dashes=False, 
        markersize=9, linewidth=2.5, palette=args.palette, ax=axes[1]
    )
    axes[1].set_title('Peak Memory Usage (Log Scale)', fontsize=15, fontweight='bold')
    axes[1].set_xlabel('Dataset (Mutation Rate / Abundance)', fontsize=13)
    axes[1].set_ylabel('Peak Memory (MB) - Log10', fontsize=13)
    axes[1].set_yscale('log') 
    axes[1].tick_params(axis='x', rotation=30)
    if axes[1].get_legend(): axes[1].get_legend().remove()

    # ==========================================
    # 🌟 全局图例处理与保存
    # ==========================================
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.02), 
               ncol=min(10, len(labels)), frameon=False, fontsize=11)

    plt.subplots_adjust(bottom=0.25, top=0.9) 
    
    out_prefix = args.out_prefix
    if out_prefix.endswith('/') or out_prefix.endswith('\\'):
        out_prefix = os.path.join(out_prefix, 'resource_evaluation')
        
    out_png = f"{out_prefix}.png"
    out_pdf = f"{out_prefix}.pdf"
    
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    
    print(f"🎉 绘图大功告成！文件已保存:\n  👉 {out_png}\n  👉 {out_pdf}")

if __name__ == "__main__":
    main()
