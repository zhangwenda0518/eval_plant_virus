#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_radar(input_dir, out_dir, metrics_list, output_format):
    os.makedirs(out_dir, exist_ok=True)
    tsv_files = glob.glob(os.path.join(input_dir, "*.tsv"))
    
    if not tsv_files:
        print(f"❌ 错误: 在 {input_dir} 目录下找不到任何 .tsv 文件！")
        return

    print(f"🚀 开始使用 Python (Plotly) 绘制 OPAL 雷达图...")
    print(f"📁 输入目录: {input_dir}")
    print(f"📁 输出目录: {out_dir}")
    print(f"📊 绘制指标: {', '.join(metrics_list)}\n")

    # Plotly 经典的高对比度色板
    colors = ['#E69F00', '#56B4E9', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#0072B2']

    success_count = 0
    for file in tsv_files:
        rank_name = os.path.splitext(os.path.basename(file))[0]
        
        try:
            df = pd.read_csv(file, sep='\t')
        except Exception as e:
            print(f"⚠️ 无法读取 {file}: {e}")
            continue
            
        # 提取当前文件中实际存在的指标
        available_metrics = [m for m in metrics_list if m in df.columns]
        if len(available_metrics) < 3:
            print(f"⏭️ 跳过 [{rank_name}]: 有效的评估指标少于 3 个，无法构成多边形雷达图。")
            continue
            
        print(f"   => 正在绘制: {rank_name.capitalize()} ...")
        
        samples = df['sample'].unique()
        tools = df['tool'].unique()
        
        # 按样本数量动态创建 1 行 N 列的分面子图
        fig = make_subplots(
            rows=1, cols=len(samples), 
            subplot_titles=samples,
            specs=[[{'type': 'polar'}] * len(samples)]
        )
        
        for i, sample in enumerate(samples):
            sample_df = df[df['sample'] == sample]
            
            for j, tool in enumerate(tools):
                tool_df = sample_df[sample_df['tool'] == tool]
                if tool_df.empty: continue
                
                # 提取数值并闭合多边形首尾相连
                values = tool_df[available_metrics].iloc[0].tolist()
                values += [values[0]] 
                theta = available_metrics + [available_metrics[0]]
                
                fig.add_trace(
                    go.Scatterpolar(
                        r=values,
                        theta=theta,
                        fill='toself',
                        name=tool,
                        line=dict(color=colors[j % len(colors)], width=2),
                        opacity=0.4,
                        showlegend=(i == 0) # 只在第一个子图显示图例
                    ),
                    row=1, col=i+1
                )
        
        # 整体排版美化
        fig.update_layout(
            title=dict(
                text=f"OPAL Performance Comparison - Rank: <b>{rank_name.capitalize()}</b>", 
                x=0.5, font=dict(size=22)
            ),
            # 锁定雷达图刻度为 0 到 1.1（留出边缘呼吸感）
            polar=dict(radialaxis=dict(visible=True, range=[0, 1.1], showticklabels=True)),
            polar2=dict(radialaxis=dict(visible=True, range=[0, 1.1], showticklabels=True)) if len(samples) > 1 else None,
            polar3=dict(radialaxis=dict(visible=True, range=[0, 1.1], showticklabels=True)) if len(samples) > 2 else None,
            polar4=dict(radialaxis=dict(visible=True, range=[0, 1.1], showticklabels=True)) if len(samples) > 3 else None,
            polar5=dict(radialaxis=dict(visible=True, range=[0, 1.1], showticklabels=True)) if len(samples) > 4 else None,
            template="plotly_white",
            legend=dict(title=dict(text="<b>Tools</b>"), orientation="v", y=0.5)
        )
        
        # 输出文件
        base_out_path = os.path.join(out_dir, f"radar_{rank_name}")
        if output_format in ['html', 'all']:
            fig.write_html(f"{base_out_path}.html")
        if output_format in ['pdf', 'all']:
            try:
                fig.write_image(f"{base_out_path}.pdf", width=1400, height=600)
            except ValueError:
                print("   [提示] 导出 PDF 需要安装 kaleido: pip install -U kaleido")
        if output_format in ['png', 'all']:
            try:
                fig.write_image(f"{base_out_path}.png", width=1400, height=600, scale=2)
            except ValueError:
                pass
                
        success_count += 1

    print(f"\n🎉 完美！成功绘制了 {success_count} 个分类等级的雷达图。")

def main():
    parser = argparse.ArgumentParser(
        description="🕸️ OPAL 多维度评估雷达图绘制工具 (Python + Plotly 版)\n自动读取 by_rank 目录下的评估结果，按样本分面绘制可交互的雷达图。",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input-dir', default="by_rank",
                        help="输入目录，包含 OPAL 评估生成的 tsv 文件 (如 class.tsv)\n(默认: by_rank)")
    parser.add_argument('-o', '--out-dir', default="python_radar_plots",
                        help="雷达图输出目录\n(默认: python_radar_plots)")
    parser.add_argument('-m', '--metrics', default="Completeness,Purity,L1 norm error,Weighted UniFrac error",
                        help="需要绘制的评估指标，用逗号分隔\n(默认: Completeness,Purity,L1 norm error,Weighted UniFrac error)")
    parser.add_argument('-f', '--format', default="html", choices=['html', 'pdf', 'png', 'all'],
                        help="输出格式。'html'为交互式网页，'pdf'/'png'为静态图 (需装kaleido)\n(默认: html)")

    args = parser.parse_args()
    
    # 将逗号分隔的字符串转换为列表
    metrics_list = [m.strip() for m in args.metrics.split(',')]
    
    plot_radar(args.input_dir, args.out_dir, metrics_list, args.format)

if __name__ == "__main__":
    main()
