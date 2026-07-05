#!/usr/bin/env python3
"""
评估三：候选病毒鉴定策略 — 完整评估脚本
支持 9 工具 × 3 过滤模式 (Raw / UniProt Filtered / Strict) × 搜索原理分组 (P1-P5) × 投票策略
新增: 3x1分组柱状图展示与雷达图多维评估
"""

import argparse, os, sys
from itertools import combinations
import pandas as pd
import numpy as np
from collections import OrderedDict
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, average_precision_score, auc

# ── 全局配置 ────────────────────────────────────────────────
PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#7E6148',
           '#BB0021', '#5C6BC0', '#FFB74D']
TOOL_NAMES = {
    'blast': 'BLASTn', 'genomad': 'geNomad', 'rdrpcatch': 'RdRpCatch',
    'viralm': 'ViraLM', 'virbot': 'VirBot', 'metabuli': 'MetaBuli',
    'viralverify': 'ViralVerify', 'virhunter': 'VirHunter', 'virsorter2': 'VirSorter2',
    'all': 'Ensemble',
}
FILTER_MODES = ['raw', 'filter', 'strict']
FILTER_LABELS = {'raw': 'Raw', 'filter': 'UniProt Filtered', 'strict': 'Strict Filtered'}

PREDEFINED_COMBOS = {
    'P1_Nucleic':     ['BLASTn'],
    'P2_Protein':     ['VirBot', 'VirSorter2', 'ViralVerify'],
    'P3_RdRp':        ['RdRpCatch'],
    'P4_DeepLearn':   ['ViraLM', 'VirHunter', 'geNomad'],
    'P5_Kmer':        ['MetaBuli'],
    'P5_Ensemble':    ['BLASTn','VirBot','VirSorter2','ViralVerify','RdRpCatch','ViraLM','VirHunter','geNomad','MetaBuli'],
}

# ── 指标计算 ────────────────────────────────────────────────

def confusion_stats(y_true, y_pred):
    return (int(np.sum((y_true==1)&(y_pred==1))), int(np.sum((y_true==0)&(y_pred==1))),
            int(np.sum((y_true==1)&(y_pred==0))), int(np.sum((y_true==0)&(y_pred==0))))

def compute_metrics(y_true, y_pred, label=''):
    TP, FP, FN, TN = confusion_stats(y_true, y_pred)
    precision = TP/max(1, TP+FP); recall = TP/max(1, TP+FN)
    f1 = 2*precision*recall/max(1e-10, precision+recall)
    accuracy = (TP+TN)/max(1, len(y_true))
    mcc_num = TP*TN - FP*FN
    mcc_den = np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    mcc = mcc_num/max(1e-10, mcc_den)
    return OrderedDict(Tool=label, N_True=int(np.sum(y_true)), N_Pred=int(np.sum(y_pred)),
                       TP=TP, FP=FP, FN=FN, TN=TN,
                       Precision=round(precision,4), Recall=round(recall,4),
                       F1=round(f1,4), Accuracy=round(accuracy,4), MCC=round(mcc,4))

def compute_auprc(y_true, y_score):
    prec, rec, _ = precision_recall_curve(y_true, y_score)
    ap = average_precision_score(y_true, y_score)
    return ap, auc(rec, prec), prec, rec

# ── 数据加载 ────────────────────────────────────────────────

def load_labels(path):
    df = pd.read_csv(path, sep='\t')
    df['Coverage'] = df['seq_id'].str.extract(r'cov(\d+)')[0].astype(float)
    print(f"[Labels] {len(df)} seqs | Pos={sum(df['label']=='positive')} | Neg={sum(df['label']!='positive')} | Types={dict(df['type'].value_counts())}")
    return df

def load_predictions(path, labels_df):
    with open(path) as f:
        predicted = set(l.strip() for l in f if l.strip())
    return labels_df['seq_id'].isin(predicted).astype(int).values

# ── 绘图功能 ────────────────────────────────────────────────

def set_style():
    sns.set_theme(style='whitegrid', font_scale=1.1)
    plt.rcParams.update({'figure.dpi':150,'savefig.dpi':300,'savefig.bbox':'tight'})

# 【新功能】3x1 分组条形图
def plot_filter_comparison(overall_df, outpath):
    if 'Mode' not in overall_df.columns: return
    set_style()
    df_clean = overall_df.copy()
    # 提取纯净工具名
    df_clean['Base_Tool'] = df_clean['Tool'].apply(lambda x: x.split(' (')[0])
    tools = sorted(df_clean['Base_Tool'].unique())
    if 'Ensemble' in tools:
        tools.remove('Ensemble'); tools.append('Ensemble') # 把集成放最后

    modes = ['Raw', 'UniProt Filtered', 'Strict Filtered']
    actual_modes = [m for m in modes if m in df_clean['Mode'].unique()]
    if not actual_modes: actual_modes = sorted(df_clean['Mode'].unique())

    metrics = ['Recall', 'Precision', 'F1']
    
    # 3行1列的布局
    fig, axes = plt.subplots(3, 1, figsize=(14, 12)) 
    x = np.arange(len(tools))
    width = 0.8 / len(actual_modes)
    colors = ['#76C8A6', '#8FAADC', '#FDE05A'] 

    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        for i, mode in enumerate(actual_modes):
            md = df_clean[df_clean['Mode'] == mode].set_index('Base_Tool')
            vals = [float(md.loc[t, metric]) if t in md.index else 0 for t in tools]
            ax.bar(x + i * width, vals, width, alpha=0.9, 
                   label=mode, color=colors[i % len(colors)], edgecolor='white', linewidth=0.5)
            
        ax.set_xticks(x + width * (len(actual_modes) - 1) / 2)
        ax.set_xticklabels(tools, rotation=20, ha='right', fontsize=11)
        ax.set_ylabel(metric, fontsize=12, fontweight='bold')
        ax.set_title(f'{metric} by Filter Mode', fontsize=14)
        ax.set_ylim(0, 1.15)
        ax.grid(alpha=0.3, axis='y', linestyle='--')
        
        if idx == 0:
            ax.legend(fontsize=11, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3)

    plt.tight_layout()
    out_file = os.path.join(outpath, 'Fig_Filter_Comparison_Combined_3x1.png')
    fig.savefig(out_file, dpi=300)
    plt.close(fig)
    print(f"  Saved: {os.path.basename(out_file)}")

# 【新功能】3种过滤模式的雷达图
def plot_radar_charts(overall_df, outpath):
    if 'Mode' not in overall_df.columns: return
    set_style()
    df_clean = overall_df.copy()
    df_clean['Base_Tool'] = df_clean['Tool'].apply(lambda x: x.split(' (')[0])
    
    modes = ['Raw', 'UniProt Filtered', 'Strict Filtered']
    actual_modes = [m for m in modes if m in df_clean['Mode'].unique()]
    metrics = ['Precision', 'Recall', 'F1']
    metric_colors = ['#4DBBD5', '#00A087', '#E64B35'] # 蓝(P), 绿(R), 红(F1)

    for mode in actual_modes:
        md = df_clean[df_clean['Mode'] == mode].set_index('Base_Tool')
        tools = sorted(md.index.tolist())
        if 'Ensemble' in tools:
            tools.remove('Ensemble'); tools.append('Ensemble')

        # 计算雷达图角度
        num_vars = len(tools)
        angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
        angles += angles[:1] # 闭合多边形

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
        
        for idx, metric in enumerate(metrics):
            values = [float(md.loc[t, metric]) if t in md.index else 0 for t in tools]
            values += values[:1] # 闭合多边形
            
            ax.plot(angles, values, color=metric_colors[idx], linewidth=2.5, label=metric)
            ax.fill(angles, values, color=metric_colors[idx], alpha=0.15)
            
        # 设置雷达图样式
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(tools, fontsize=11)
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], color="grey", size=8)
        
        ax.set_title(f'Performance Radar\n[{mode}]', size=15, fontweight='bold', y=1.1)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
        
        out_file = os.path.join(outpath, f'Fig_Radar_{mode.replace(" ", "_")}.png')
        plt.tight_layout()
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        print(f"  Saved: {os.path.basename(out_file)}")

def plot_pr_curves(pr_data, outpath):
    set_style(); fig, ax = plt.subplots(figsize=(8,7))
    for i,(label,(prec,rec,ap)) in enumerate(pr_data.items()):
        ax.plot(rec,prec,linewidth=2,color=PALETTE[i%len(PALETTE)],label=f'{label} (AUPRC={ap:.3f})')
    ax.set_xlabel('Recall'); ax.set_ylabel('Precision'); ax.set_title('PR Curves',fontweight='bold')
    ax.legend(loc='lower left',fontsize=8); ax.set_xlim(0,1.05); ax.set_ylim(0,1.05); ax.grid(alpha=0.3)
    fig.savefig(outpath); plt.close(fig); print(f"  Saved: {os.path.basename(outpath)}")

def plot_fp_reduction(overall_df, outpath):
    if 'Mode' not in overall_df.columns or 'FP' not in overall_df.columns: return
    set_style(); tools = sorted(overall_df['Tool'].apply(lambda x: x.split(' (')[0]).unique())
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,6))
    
    for i, tool in enumerate(tools):
        sub = overall_df[overall_df['Tool'].str.startswith(tool)]
        vals, fp_rates, modes_found = [], [], []
        
        raw_fp = sub[sub['Mode']==FILTER_LABELS['raw']]['FP'].values
        raw_fp = int(raw_fp[0]) if len(raw_fp) else 1
        
        for m in FILTER_MODES:
            row = sub[sub['Mode']==FILTER_LABELS[m]]
            if len(row):
                fp = int(row['FP'].values[0])
                vals.append(fp); fp_rates.append(fp/max(1,raw_fp)*100)
                modes_found.append(FILTER_LABELS[m])
                
        if vals:
            ax1.plot(modes_found, vals, marker='o', linewidth=2, color=PALETTE[i%len(PALETTE)], label=tool)
            ax2.plot(modes_found, fp_rates, marker='s', linewidth=2, color=PALETTE[i%len(PALETTE)], label=tool)
            
    ax1.set_ylabel('False Positives'); ax1.set_title('FP Count by Filter Mode')
    ax1.legend(fontsize=8, bbox_to_anchor=(1.05,1)); ax1.grid(alpha=0.3)
    ax2.set_ylabel('FP Relative to Raw (%)'); ax2.set_title('FP Reduction Rate')
    ax2.legend(fontsize=8, bbox_to_anchor=(1.05,1)); ax2.grid(alpha=0.3); ax2.set_ylim(0,120)
    
    plt.tight_layout(); fig.savefig(os.path.join(outpath,'Fig_FP_Reduction.png'),dpi=300)
    plt.close(fig); print(f"  Saved: Fig_FP_Reduction.png")

def plot_overlap(y_preds_dict, outpath):
    set_style(); fig, axes = plt.subplots(1,2,figsize=(14,6))
    tools_bar = [t for t in y_preds_dict.keys() if t!='Ensemble']
    unique = [int(np.sum(y_preds_dict[t] & ~np.logical_or.reduce([y_preds_dict[s] for s in tools_bar if s!=t]))) for t in tools_bar]
    shared = [int(np.sum(y_preds_dict[t]))-u for t,u in zip(tools_bar,unique)]
    x = np.arange(len(tools_bar)); width = 0.35
    axes[0].bar(x-width/2,unique,width,label='Unique',color='#E64B35',alpha=0.85)
    axes[0].bar(x+width/2,shared,width,label='Shared',color='#4DBBD5',alpha=0.85)
    axes[0].set_xticks(x); axes[0].set_xticklabels(tools_bar,rotation=45,ha='right')
    axes[0].set_ylabel('Predicted'); axes[0].set_title('Overlap: Unique vs Shared')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3,axis='y')
    votes = np.sum(np.column_stack([y_preds_dict[t] for t in tools_bar]),axis=1)
    n_positive = len(tools_bar)
    counts = [np.sum(votes==i) for i in range(1,n_positive+1)]
    axes[1].bar(range(1,n_positive+1),counts,color=plt.cm.YlOrRd(np.linspace(0.3,0.9,n_positive)),edgecolor='white')
    axes[1].set_xlabel('Tools detecting'); axes[1].set_ylabel('Sequences')
    axes[1].set_title('Detection Consensus'); axes[1].grid(alpha=0.3,axis='y')
    for i,v in enumerate(counts):
        axes[1].text(i+1,v+1,str(v),ha='center',fontsize=10)
    plt.tight_layout(); fig.savefig(os.path.join(outpath,'Fig_Detection_Overlap.png'),dpi=300)
    plt.close(fig); print(f"  Saved: Fig_Detection_Overlap.png")

# ── 核心逻辑 ───────────────────────────────────────────────────

def process_one_mode(labels, args, mode):
    mode_label = FILTER_LABELS[mode]
    if mode == 'raw':
        sub_dir = args.result_dir; suffix = args.suffix
    elif mode == 'filter':
        sub_dir = os.path.join(args.result_dir, 'uniprot_filter_output_filter'); suffix = '.uniprot_filtered.id'
    else: 
        sub_dir = os.path.join(args.result_dir, 'uniprot_filter_output_strict'); suffix = '.uniprot_filtered.id'

    tool_keys = args.tools.split(',')
    if args.ensemble: tool_keys.append('all')

    all_overall = []; ds_store = {}; y_preds_out = {}

    for tool in tool_keys:
        fname = f'{args.prefix}.{tool}{suffix}' if tool != 'all' else f'{args.prefix}.all{suffix}'
        fpath = os.path.join(sub_dir, fname)
        if not os.path.exists(fpath): continue

        base_label = TOOL_NAMES.get(tool, tool.upper())
        tool_label = base_label if mode == 'raw' else f"{base_label} ({mode_label})"
        
        metrics, y_pred = evaluate_tool(fpath, labels, tool_label)
        metrics['Mode'] = mode_label; all_overall.append(metrics)
        
        y_true = (labels['label']=='positive').astype(int).values
        ap, _, prec, rec = compute_auprc(y_true, y_pred.astype(float))
        ds_store[tool_label] = (prec, rec, ap)
        all_overall[-1]['AUPRC'] = round(ap, 4)
        y_preds_out[tool_label] = y_pred

    overall_df = pd.DataFrame(all_overall)
    overall_no_cov = overall_df[overall_df['Precision'].notna()].drop_duplicates(subset='Tool') if not overall_df.empty else pd.DataFrame()
    return overall_no_cov, ds_store, y_preds_out

def main():
    parser = argparse.ArgumentParser(description='病毒鉴定评估脚本 (带多维图表)')
    parser.add_argument('--result-dir', default='step8_result/step3_identification_eval')
    parser.add_argument('--labels', required=True)
    parser.add_argument('--outdir', default='step8_result_analysis')
    parser.add_argument('--tools', default='blast,genomad,rdrpcatch,viralm,virbot,metabuli,viralverify,virhunter,virsorter2')
    parser.add_argument('--ensemble', action='store_true', default=True)
    parser.add_argument('--prefix', default='step3_identification_eval_virus')
    parser.add_argument('--suffix', default='.result.id')
    parser.add_argument('--filter-mode', default='raw', choices=['raw','filter','strict','all'])
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    labels = load_labels(args.labels)
    modes_to_run = FILTER_MODES if args.filter_mode == 'all' else [args.filter_mode]

    all_overall_across = []; ds_stores = {}; y_preds_per_mode = {}

    for mode in modes_to_run:
        overall, ds_store, y_preds_out = process_one_mode(labels, args, mode)
        if not overall.empty:
            all_overall_across.append(overall)
            ds_stores[mode] = ds_store
            y_preds_per_mode[mode] = y_preds_out
            overall.to_csv(os.path.join(args.outdir, f'identification_overall_{mode}.tsv'), sep='\t', index=False)

    if not all_overall_across:
        print("未找到任何结果文件，退出。"); return

    raw_overall = all_overall_across[modes_to_run.index('raw')] if 'raw' in modes_to_run else all_overall_across[0]

    # ── 制图阶段 ──
    if not args.no_plot:
        # 1. 绘制综合表格 (包含 Accuracy 和 MCC)
        tbl_cols = [c for c in ['Tool','Precision','Recall','F1','AUPRC','MCC','Accuracy'] if c in raw_overall.columns]
        fig_tbl, ax_tbl = plt.subplots(figsize=(10, 3)); ax_tbl.axis('off')
        tbl_data = raw_overall[tbl_cols].round(4)
        tbl = ax_tbl.table(cellText=tbl_data.values, colLabels=tbl_data.columns, cellLoc='center', loc='center')
        tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1.3, 1.8)
        fig_tbl.savefig(os.path.join(args.outdir, 'Fig_Metrics_Summary_Table.png'), dpi=300, bbox_inches='tight')
        plt.close(fig_tbl); print(f"  Saved: Fig_Metrics_Summary_Table.png")

        # 2. PR 曲线与重叠图
        raw_ds = ds_stores.get('raw', ds_stores.get(FILTER_MODES[0], {}))
        if raw_ds:
            plot_pr_curves(raw_ds, os.path.join(args.outdir, 'Fig_Identification_PR_Curves.png'))
            raw_yp = y_preds_per_mode.get('raw', {})
            if raw_yp: plot_overlap(raw_yp, args.outdir)

        # 3. 多模式对比图 (仅在全模式时绘制)
        if len(modes_to_run) > 1:
            combined = pd.concat(all_overall_across)
            combined.to_csv(os.path.join(args.outdir, 'identification_filter_comparison.tsv'), sep='\t', index=False)
            plot_fp_reduction(combined, args.outdir)
            
            # 【应用新布局】3x1 分组条形图 (Recall, Precision, F1)
            plot_filter_comparison(combined, args.outdir)
            
            # 【应用雷达图】按模式分开生成
            plot_radar_charts(combined, args.outdir)

    print(f"\n🎉 评估完成: {args.outdir}")

if __name__ == '__main__':
    main()
