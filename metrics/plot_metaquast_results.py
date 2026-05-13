#!/usr/bin/env python3
"""MetaQUAST 7组结果汇总绘图"""
import argparse, os, glob
import pandas as pd, numpy as np
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

PALETTE = ['#4C72B0','#DD8452','#55A868','#C44E52','#8B5CF6','#F39B7F','#95A5A6']
TOOL_ORDER = ['Megahit','RNAViralSPAdes','Penguin',
              'MH_Merge','MH_SplitMerge','ALL_Merge','ALL_SplitMerge']

def parse_all(dir_path):
    """解析所有 transposed_report.tsv"""
    pattern = os.path.join(dir_path, "*", "runs_per_reference", "*", "transposed_report.tsv")
    dfs = []
    for f in glob.glob(pattern):
        try:
            df = pd.read_csv(f, sep='\t')
            parts = f.split(os.sep)
            df['Sample'] = parts[-3]
            df['Virus'] = parts[-2]
            df['Tool'] = df['Assembly'].str.replace('_contig','').str.replace('_refineC','')
            dfs.append(df)
        except: continue
    return pd.concat(dfs) if dfs else pd.DataFrame()

def plot_quality(quality_df, outdir):
    """基因组覆盖度+NGA50 箱线图"""
    if quality_df.empty: return
    q = quality_df.copy()
    q['Tool'] = q['Tool'].replace({'MH_merge':'MH_Merge','MH_split_merge':'MH_SplitMerge',
                                     'ALL_merge':'ALL_Merge','all_tools_refineC_merge':'ALL_SplitMerge'})
    q = q[q['Tool'].isin(TOOL_ORDER)]

    metrics = [('Genome fraction (%)','Fig_MQ_GenomeFraction','Genome Fraction (%)'),
               ('NGA50','Fig_MQ_NGA50','NGA50 (bp)')]
    for col, fname, ylabel in metrics:
        if col not in q.columns: continue
        fig, ax = plt.subplots(figsize=(10,6))
        sns.boxplot(data=q, x='Tool', y=col, order=TOOL_ORDER, palette=PALETTE, ax=ax)
        ax.set_ylabel(ylabel, fontsize=13)
        ax.set_title(f'Assembly Quality: {ylabel}', fontsize=14, fontweight='bold')
        ax.tick_params(axis='x', rotation=20)
        ax.grid(axis='y', alpha=0.3)
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f'{fname}.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}.png")

def plot_resource_summary(asm_dir, log_dir, outdir):
    """7组资源消耗汇总柱状图"""
    # 从 time.mem.log 提取数据
    pass  # 简略实现，需要时补充

def main():
    parser = argparse.ArgumentParser('MetaQUAST 结果绘图')
    parser.add_argument('--metaquast-dir', required=True, help='MetaQUAST输出目录')
    parser.add_argument('--asm-dir', help='组装输出目录(用于资源图)')
    parser.add_argument('--log-dir', help='日志目录(用于资源图)')
    parser.add_argument('--outdir', required=True, help='图片输出目录')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    quality_df = parse_all(args.metaquast_dir)
    if not quality_df.empty:
        quality_df.to_csv(os.path.join(args.outdir, 'metaquast_quality.tsv'), sep='\t', index=False)
        plot_quality(quality_df, args.outdir)
    print(f"Done. Output: {args.outdir}")

if __name__ == '__main__':
    main()
