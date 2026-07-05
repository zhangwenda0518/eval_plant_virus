#!/usr/bin/env python3
"""
评估三：候选病毒鉴定策略 — 完整且全功能版 (Polars加速 + 实体级评估 + FP来源细化)
支持 9 工具 × 3 过滤模式 × 搜索原理分组 (P1-P5) × 投票策略
新增: 多片段病毒实体评估、假阳性来源堆叠图、3x1多模式柱状图、雷达图

用法:
  python eval_identification.py --result-dir step8_result/step3_identification_eval \
      --labels step3_identification_eval/sequence_labels.tsv --outdir analysis/ --filter-mode all

输出目录结构:
  {outdir}/
    overall/      — 总体评估: PR曲线、ABCD四合一、汇总表、雷达图 + TSV
    modes/        — 多模式对比: FP减少、Precision-Recall散点、3x1柱状图 + TSV
    overlap/      — 工具重叠分析 + TSV
    combinations/ — 组合分析: 热力图、Top20 + TSV
    coverage/     — 覆盖率召回曲线 + TSV
    entity/       — 实体级(分段病毒)评估 + 假阳性来源堆叠图 + TSV
"""

import argparse, os, re, sys
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from collections import OrderedDict
from Bio import SeqIO
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, average_precision_score, auc

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

# ── 全局配置 ────────────────────────────────────────────────
PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
           '#7E6148', '#BB0021', '#5C6BC0', '#FFB74D', '#8C564B']
TOOL_NAMES = {
    'blast': 'BLASTn', 'genomad': 'geNomad', 'rdrpcatch': 'RdRpCatch',
    'viralm': 'ViraLM', 'virbot': 'VirBot', 'metabuli': 'MetaBuli',
    'viralverify': 'ViralVerify', 'virhunter': 'VirHunter', 'virsorter2': 'VirSorter2',
    'all': 'Ensemble',
}
# 全局工具-颜色锚定: 确保同一工具在所有图中颜色绝对一致
_ALL_TOOLS_LIST = sorted([v for k, v in TOOL_NAMES.items() if k != 'all'])
_ALL_TOOLS_LIST.append('Ensemble')
TOOL_COLORS = {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(_ALL_TOOLS_LIST)}

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

# ── Polars 加速数据加载 ─────────────────────────────────────

class DataEngine:
    """
    统一数据引擎: 优先使用 Polars 加速, 降级到 Pandas
    对外暴露一致的 DataFrame-like API
    """
    def __init__(self):
        self.use_polars = HAS_POLARS
        self._df = None
        if self.use_polars:
            print("[Engine] 启用 Polars 加速模式")
        else:
            print("[Engine] Polars 未安装, 降级到 Pandas")

    def load_labels(self, path):
        if self.use_polars:
            df = pl.read_csv(path, separator='\t')
            df = df.with_columns(
                pl.col("seq_id").str.extract(r'cov(\d+)').cast(pl.Float64).alias("Coverage")
            )
            if "Category" not in df.columns:
                df = df.with_columns(pl.col("type").alias("Category"))
            self._df = df
            pos = df.filter(pl.col("label") == "positive").height
            neg = df.filter(pl.col("label") != "positive").height
            print(f"[Labels] {df.height} seqs | Pos={pos} | Neg={neg}")
            return df
        else:
            df = pd.read_csv(path, sep='\t')
            df['Coverage'] = df['seq_id'].str.extract(r'cov(\d+)').astype(float)
            if 'Category' not in df.columns:
                df['Category'] = df['type']
            self._df = df
            pos = sum(df['label'] == 'positive')
            neg = sum(df['label'] != 'positive')
            print(f"[Labels] {len(df)} seqs | Pos={pos} | Neg={neg}")
            return df

    def load_predictions(self, path):
        with open(path) as f:
            predicted = set(l.strip() for l in f if l.strip())
        if self.use_polars:
            return self._df.select(pl.col("seq_id").is_in(predicted).cast(pl.Int32)).to_series().to_numpy()
        else:
            return self._df['seq_id'].isin(predicted).astype(int).values

    def get_y_true(self):
        if self.use_polars:
            return (self._df.get_column("label") == "positive").cast(pl.Int32).to_numpy()
        else:
            return (self._df['label'] == 'positive').astype(int).values

    def filter_short_source_viruses(self, virus_dir, min_length):
        """
        排除源病毒基因组短于 min_length bp 的阳性序列。
        类病毒(<1000 bp)的鉴定本质上是不同的问题, 应单独评估。
        """
        if not virus_dir or not os.path.isdir(virus_dir) or min_length <= 0:
            return self._df

        # 1. 扫描病毒fasta获取长度
        virus_lengths = {}
        for f in Path(virus_dir).glob("*.fasta"):
            try:
                for rec in SeqIO.parse(f, "fasta"):
                    acc = rec.id.split()[0]
                    virus_lengths[acc] = len(rec.seq)
            except Exception:
                continue

        short_sources = {acc for acc, length in virus_lengths.items() if length < min_length}
        if not short_sources:
            print(f"[filter] 所有病毒 ≥ {min_length} bp，无需过滤")
            return self._df

        # 2. 提取 seq_id 中的 source accession
        if self.use_polars:
            df = self._df.with_columns(
                pl.col("seq_id").str.extract(r'source=([^|]+)').alias("_source")
            )
            before = df.height
            df = df.filter(
                (pl.col("label") != "positive") |
                (~pl.col("_source").is_in(short_sources))
            )
            df = df.drop("_source")
            after = df.height
            self._df = df
        else:
            df = self._df.copy()
            df['_source'] = df['seq_id'].str.extract(r'source=([^|]+)')
            before = len(df)
            df = df[~((df['label'] == 'positive') & (df['_source'].isin(short_sources)))]
            df = df.drop(columns=['_source'])
            after = len(df)
            self._df = df

        pos_removed = before - after
        pos_after = self._df.filter(pl.col("label") == "positive").height if self.use_polars else sum(self._df['label'] == 'positive')
        neg_after = self._df.filter(pl.col("label") != "positive").height if self.use_polars else sum(self._df['label'] != 'positive')
        print(f"[filter] 排除 {len(short_sources)} 个短源病毒 (< {min_length} bp): {sorted(short_sources)}")
        print(f"[filter] 移除 {pos_removed} 条阳性序列")
        print(f"[Labels] {after} seqs | Pos={pos_after} | Neg={neg_after}")
        return self._df

    @property
    def df(self):
        """返回原始 DataFrame (Polars 或 Pandas)"""
        return self._df


# ── 切片评估函数 (与引擎解耦) ────────────────────────────────

def evaluate_tool(pred_path, labels_df, labels_pd, engine, tool_label=''):
    """评估单个工具的预测结果 (兼容两种后端)"""
    y_pred = engine.load_predictions(pred_path)
    y_true = engine.get_y_true()
    return compute_metrics(y_true, y_pred, tool_label), y_pred

def by_type_eval(labels_pd, y_pred, tool_label):
    """按序列类型分层评估"""
    rows = []
    for stype in sorted(labels_pd['type'].unique()):
        mask = labels_pd['type'] == stype
        if mask.sum() == 0: continue
        yt = (labels_pd.loc[mask, 'label'] == 'positive').astype(int).values
        yp = y_pred[mask.values]
        rows.append({'Type': stype, **compute_metrics(yt, yp, tool_label)})
    return rows

def by_coverage_eval(labels_pd, y_pred, tool_label):
    """按覆盖率分层评估正样本召回率"""
    rows = []
    pos = labels_pd[labels_pd['label'] == 'positive']
    for cov in sorted(pos['Coverage'].dropna().unique()):
        mask = pos['Coverage'] == cov
        n = int(mask.sum())
        detected = int(np.sum(y_pred[pos[mask].index]))
        rows.append({
            'Coverage': f'{int(cov)}%', 'N_Seqs': n, 'Detected': detected,
            'Recall': round(detected / max(1, n), 4), 'Tool': tool_label
        })
    return rows

def _calc_metrics(tp, fp, fn, tn):
    """从混淆矩阵计算 P/R/F1/MCC/Accuracy"""
    p = tp / max(1, tp + fp)
    r = tp / max(1, tp + fn)
    f1 = 2 * p * r / max(1e-10, p + r)
    mcc_num = tp * tn - fp * fn
    mcc_den = np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = mcc_num / max(1e-10, mcc_den)
    acc = (tp + tn) / max(1, tp + fp + fn + tn)
    return round(p,4), round(r,4), round(f1,4), round(mcc,4), round(acc,4)

def _calc_auprc_stratum(pos_mask_full, y_pred, y_true_all):
    """计算特定正样本层的AUPRC。pos_mask_full: 与y_true_all等长的布尔数组"""
    yt = np.zeros(len(y_true_all), dtype=int)
    mask_arr = pos_mask_full if isinstance(pos_mask_full, np.ndarray) else pos_mask_full.values
    yt[mask_arr] = 1
    yp = y_pred.astype(float)
    ap, _, _, _ = compute_auprc(yt, yp)
    return round(ap, 4)

def by_scheme_stratified_eval(labels_pd, y_pred, tool_label):
    """
    MASTER数据集分层评估: scheme (A/B/C) × mutation × coverage
    计算 Precision / Recall / F1 / MCC / AUPRC (包含全量负样本)
    """
    pos = labels_pd[labels_pd['label'] == 'positive'].copy()
    pos['pred'] = y_pred[pos.index]
    neg_mask = labels_pd['label'] != 'positive'
    neg_pred = y_pred[neg_mask.values]
    n_neg = int(neg_mask.sum())
    y_true_all = (labels_pd['label'] == 'positive').astype(int).values

    pos['scheme'] = pos['seq_id'].str.extract(r'scheme_([A-C])')
    pos['mut'] = pos['seq_id'].str.extract(r'mut(\d+)').astype(float)
    pos['cov'] = pos['seq_id'].str.extract(r'cov(\d+)').astype(float)

    rows = []
    for scheme in sorted(pos['scheme'].dropna().unique()):
        mask_s = pos['scheme'] == scheme
        tp = int(pos.loc[mask_s, 'pred'].sum())
        fn = int(mask_s.sum()) - tp
        fp = int(neg_pred.sum()); tn = n_neg - fp
        p, r, f1, mcc, acc = _calc_metrics(tp, fp, fn, tn)
        # 构建与 labels_pd 等长的全量mask
        full_mask = labels_pd['seq_id'].isin(pos.loc[mask_s, 'seq_id']).values
        auprc = _calc_auprc_stratum(full_mask, y_pred, y_true_all)
        rows.append({
            'Tool': tool_label, 'Scheme': f'Scheme_{scheme}', 'Stratum': 'Overall',
            'N_Pos': int(mask_s.sum()), 'N_Neg': n_neg, 'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
            'Precision': p, 'Recall': r, 'F1': f1, 'MCC': mcc, 'Accuracy': acc, 'AUPRC': auprc,
        })

        if scheme == 'A':
            for mut in sorted(pos.loc[mask_s, 'mut'].dropna().unique()):
                mask_m = mask_s & (pos['mut'] == mut)
                tp_m = int(pos.loc[mask_m, 'pred'].sum())
                fn_m = int(mask_m.sum()) - tp_m
                p_m, r_m, f1_m, mcc_m, acc_m = _calc_metrics(tp_m, fp, fn_m, tn)
                rows.append({
                    'Tool': tool_label, 'Scheme': 'Scheme_A', 'Stratum': f'mut_{int(mut)}%',
                    'N_Pos': int(mask_m.sum()), 'N_Neg': n_neg, 'TP': tp_m, 'FP': fp, 'FN': fn_m, 'TN': tn,
                    'Precision': p_m, 'Recall': r_m, 'F1': f1_m, 'MCC': mcc_m, 'Accuracy': acc_m, 'AUPRC': 0.0,
                })
        elif scheme == 'B':
            for cov in sorted(pos.loc[mask_s, 'cov'].dropna().unique()):
                mask_c = mask_s & (pos['cov'] == cov)
                tp_c = int(pos.loc[mask_c, 'pred'].sum())
                fn_c = int(mask_c.sum()) - tp_c
                p_c, r_c, f1_c, mcc_c, acc_c = _calc_metrics(tp_c, fp, fn_c, tn)
                rows.append({
                    'Tool': tool_label, 'Scheme': 'Scheme_B', 'Stratum': f'cov_{int(cov)}',
                    'N_Pos': int(mask_c.sum()), 'N_Neg': n_neg, 'TP': tp_c, 'FP': fp, 'FN': fn_c, 'TN': tn,
                    'Precision': p_c, 'Recall': r_c, 'F1': f1_c, 'MCC': mcc_c, 'Accuracy': acc_c, 'AUPRC': 0.0,
                })
        elif scheme == 'C':
            scheme_c = pos[mask_s]
            reps = scheme_c['seq_id'].str.extract(r'rep(\d+)')[0].dropna().unique()
            for rep in sorted(reps):
                mask_r = scheme_c['seq_id'].str.contains(f'rep{rep}')
                tp_r = int(scheme_c.loc[mask_r, 'pred'].sum())
                fn_r = int(mask_r.sum()) - tp_r
                p_r, r_r, f1_r, mcc_r, acc_r = _calc_metrics(tp_r, fp, fn_r, tn)
                rows.append({
                    'Tool': tool_label, 'Scheme': 'Scheme_C', 'Stratum': f'rep_{rep}',
                    'N_Pos': int(mask_r.sum()), 'N_Neg': n_neg, 'TP': tp_r, 'FP': fp, 'FN': fn_r, 'TN': tn,
                    'Precision': p_r, 'Recall': r_r, 'F1': f1_r, 'MCC': mcc_r, 'Accuracy': acc_r, 'AUPRC': 0.0,
                })
    return rows


def plot_scheme_stratified(strat_df, outpath, metric='Recall'):
    """3×1分层: (A)突变→指标 (B)覆盖率→指标 (C)VIROMOCK均值±std, metric∈{Recall,Precision,F1}"""
    if strat_df.empty or metric not in strat_df.columns: return
    set_style(0.9)
    fig, axes = plt.subplots(3, 1, figsize=(14, 18))
    m = metric

    # AUPRC仅Overall层有值, 改用柱状图
    if metric == 'AUPRC':
        overall = strat_df[strat_df['Stratum'] == 'Overall']
        for idx, scheme in enumerate(['Scheme_A', 'Scheme_B', 'Scheme_C']):
            ax = axes[idx]
            sub = overall[overall['Scheme'] == scheme]
            if sub.empty: continue
            tools = sorted(sub['Tool'].unique())
            vals = [float(sub[sub['Tool']==t][m].values[0]) if t in sub['Tool'].values else 0 for t in tools]
            colors = [TOOL_COLORS.get(t.split(' (')[0], '#333') for t in tools]
            ax.bar(range(len(tools)), vals, color=colors, edgecolor='white')
            ax.set_xticks(range(len(tools))); ax.set_xticklabels(tools, rotation=45, ha='right', fontsize=10)
            ax.set_ylabel(m, fontsize=13)
            labels = {'Scheme_A': 'A: Known Virus', 'Scheme_B': 'B: Novel Virus', 'Scheme_C': 'C: VIROMOCK'}
            ax.set_title(labels.get(scheme, scheme), fontsize=13, fontweight='bold')
            ax.set_ylim(0, 1.05); ax.grid(alpha=0.3, axis='y', linestyle='--')
        fig.suptitle(f'{m} by Simulation Scheme', fontsize=17, fontweight='bold', y=0.99)
        fig.subplots_adjust(top=0.95, hspace=0.4)
        fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
        print(f"  [overall] {os.path.basename(outpath)}")
        return

    # (A) Scheme A: mutation rate → recall
    ax1 = axes[0]
    mut_data = strat_df[strat_df['Stratum'].str.startswith('mut_') & ~strat_df['Stratum'].str.contains('cov')]
    if len(mut_data):
        for tool in sorted(mut_data['Tool'].unique()):
            sub = mut_data[mut_data['Tool'] == tool].copy()
            sub['mut_val'] = sub['Stratum'].str.extract(r'mut_(\d+)').astype(int)
            sub = sub.sort_values('mut_val')
            ax1.plot(sub['mut_val'], sub[m].astype(float), marker='o', linewidth=2.5,
                     color=TOOL_COLORS.get(tool.split(' (')[0], '#333'), label=tool, markersize=8)
        ax1.set_xlabel('Mutation Rate (%)', fontsize=13)
        ax1.set_ylabel(m, fontsize=13)
        ax1.set_title('(A) Scheme A: Known Virus Mutation Tolerance', fontsize=15, fontweight='bold')
        ax1.set_ylim(0, 1.05); ax1.grid(alpha=0.3, linestyle='--')
        ax1.legend(fontsize=9, ncol=5, loc='lower center', bbox_to_anchor=(0.5, -0.25), frameon=True)

    # (B) Scheme B: coverage → recall
    ax2 = axes[1]
    cov_data = strat_df[strat_df['Stratum'].str.startswith('cov_')]
    if len(cov_data):
        for tool in sorted(cov_data['Tool'].unique()):
            sub = cov_data[cov_data['Tool'] == tool].copy()
            sub['cov_val'] = sub['Stratum'].str.extract(r'cov_(\d+)').astype(int)
            sub = sub.sort_values('cov_val')
            ax2.plot(sub['cov_val'], sub[m].astype(float), marker='s', linewidth=2.5,
                     color=TOOL_COLORS.get(tool.split(' (')[0], '#333'), label=tool, markersize=8)
        ax2.set_xlabel('Coverage Level (%)', fontsize=13)
        ax2.set_ylabel(m, fontsize=13)
        ax2.set_title('(B) Scheme B: Novel Virus Coverage Sensitivity', fontsize=15, fontweight='bold')
        ax2.set_ylim(0, 1.05); ax2.grid(alpha=0.3, linestyle='--')
        ax2.legend(fontsize=9, ncol=5, loc='lower center', bbox_to_anchor=(0.5, -0.25), frameon=True)

    # (C) Scheme C: VIROMOCK — 均值折线 + std阴影
    ax3 = axes[2]
    rep_data = strat_df[strat_df['Stratum'].str.startswith('rep_')]
    if len(rep_data):
        # X轴用假序号 (rep1-5), 对每个工具计算均值±std
        tools_sorted = sorted(rep_data['Tool'].unique())
        rep_nums = sorted(rep_data['Stratum'].str.extract(r'rep_(\d+)')[0].dropna().astype(int).unique())
        for tool in tools_sorted:
            sub = rep_data[rep_data['Tool'] == tool].copy()
            sub['rep_num'] = sub['Stratum'].str.extract(r'rep_(\d+)').astype(int)
            sub = sub.sort_values('rep_num')
            color = TOOL_COLORS.get(tool.split(' (')[0], '#333')
            means, stds = [], []
            for r in rep_nums:
                vals = sub[sub['rep_num'] == r][m].astype(float)
                means.append(vals.mean() if len(vals) else 0)
                stds.append(vals.std() if len(vals) else 0)
            ax3.plot(rep_nums, means, marker='D', linewidth=2.5, color=color, label=tool, markersize=8)
            ax3.fill_between(rep_nums, [m-s for m,s in zip(means,stds)],
                             [m+s for m,s in zip(means,stds)], color=color, alpha=0.12)
        ax3.set_xlabel('VIROMOCK Replicate', fontsize=13)
        ax3.set_ylabel(f'{m} (mean ± std)', fontsize=13)
        ax3.set_title('(C) Scheme C: VIROMOCK Consistency (70% nt similarity)', fontsize=15, fontweight='bold')
        ax3.set_xticks(rep_nums)
        ax3.set_ylim(0, 1.05); ax3.grid(alpha=0.3, linestyle='--')
        ax3.legend(fontsize=9, ncol=5, loc='lower center', bbox_to_anchor=(0.5, -0.25), frameon=True)

    fig.suptitle(f'Stratified {m} by Simulation Scheme', fontsize=17, fontweight='bold', y=0.99)
    fig.subplots_adjust(top=0.95, hspace=0.4)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(outpath)}")


def plot_scheme_radar_3x1(strat_df, outpath):
    """3×1雷达图: Scheme A/B/C 三列, P/R/F1/MCC/AUPRC 五轴"""
    overall = strat_df[strat_df['Stratum'] == 'Overall']
    if overall.empty: return
    set_style(0.8)
    schemes = ['Scheme_A', 'Scheme_B', 'Scheme_C']
    scheme_labels = {'Scheme_A': 'A: Known Virus', 'Scheme_B': 'B: Novel Virus', 'Scheme_C': 'C: VIROMOCK'}
    metrics = ['Precision', 'Recall', 'F1', 'MCC', 'AUPRC']
    metric_colors = {'Precision': '#4DBBD5', 'Recall': '#00A087', 'F1': '#E64B35',
                     'MCC': '#F39B7F', 'AUPRC': '#9B59B6'}
    present_schemes = [s for s in schemes if s in overall['Scheme'].values]

    fig, axes = plt.subplots(1, len(present_schemes), figsize=(8*len(present_schemes), 7.5),
                              subplot_kw=dict(polar=True))
    if len(present_schemes) == 1: axes = [axes]

    for ax, scheme in zip(axes, present_schemes):
        sub = overall[overall['Scheme'] == scheme]
        tools = sorted(sub['Tool'].unique())
        angles = np.linspace(0, 2*np.pi, len(tools), endpoint=False).tolist()
        angles += angles[:1]

        for metric in metrics:
            values = [float(sub[sub['Tool']==t][metric].values[0]) if t in sub['Tool'].values else 0 for t in tools]
            values += values[:1]
            c = metric_colors[metric]
            ax.plot(angles, values, color=c, linewidth=2, label=metric)
            ax.fill(angles, values, color=c, alpha=0.08)
            for a, v in zip(angles[:-1], values[:-1]):
                ax.annotate(f'{v:.2f}', xy=(a, v), fontsize=7, ha='center', va='bottom', color=c, fontweight='bold')

        ax.set_xticks(angles[:-1]); ax.set_xticklabels(tools, fontsize=10)
        ax.tick_params(axis='x', pad=15)
        ax.set_ylim(0, 1.05); ax.set_yticks([0.2,0.4,0.6,0.8,1.0])
        ax.set_yticklabels(['0.2','0.4','0.6','0.8','1.0'], color='grey', size=8)
        ax.set_title(scheme_labels.get(scheme, scheme), fontsize=14, fontweight='bold', pad=25)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=9, framealpha=0.9)

    fig.suptitle('Performance Comparison Across Simulation Schemes', fontsize=15, fontweight='bold', y=1.02)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(outpath)}")


# ── 新增: 实体级评估与假阳性细化 ─────────────────────────────

def entity_level_eval(labels_pd, y_pred, tool_label):
    """
    【新增】基于 Category 的多片段病毒实体级评估
    一物种的多条序列 = 一个生物实体, 评估工具是否能完整检出所有节段
    """
    pos = labels_pd[labels_pd['label'] == 'positive'].copy()
    pos['pred'] = y_pred[pos.index]

    grouped = pos.groupby('Category').agg(
        total_segments=('pred', 'count'),
        detected_segments=('pred', 'sum')
    ).reset_index()

    total_entities = len(grouped)
    strict_hits = int(np.sum(grouped['total_segments'] == grouped['detected_segments']))
    lax_hits = int(np.sum(grouped['detected_segments'] > 0))

    # 额外: 按节段数分层
    by_nsegs = []
    for nseg_bin, label in [(1, '1'), (2, '2-3'), (4, '4+')]:
        if label == '1':
            subset = grouped[grouped['total_segments'] == 1]
        elif label == '2-3':
            subset = grouped[(grouped['total_segments'] >= 2) & (grouped['total_segments'] <= 3)]
        else:
            subset = grouped[grouped['total_segments'] >= 4]
        if len(subset) == 0:
            continue
        strict = int(np.sum(subset['total_segments'] == subset['detected_segments']))
        lax = int(np.sum(subset['detected_segments'] > 0))
        by_nsegs.append({
            'Tool': tool_label, 'Segment_Bin': label,
            'N_Entities': len(subset),
            'Strict_Recall': round(strict / max(1, len(subset)), 4),
            'Lax_Recall': round(lax / max(1, len(subset)), 4),
        })

    return {
        'Tool': tool_label,
        'Total_Entities': total_entities,
        'Complete_Entities': strict_hits,
        'Partial_Entities': lax_hits - strict_hits,
        'Missed_Entities': total_entities - lax_hits,
        'Strict_Entity_Recall': round(strict_hits / max(1, total_entities), 4),
        'Lax_Entity_Recall': round(lax_hits / max(1, total_entities), 4),
    }, by_nsegs

# 负样本子类型 → 可读标签映射
# 负样本类型标签映射 (兼容新旧两种格式)
NEG_TYPE_LABELS = {
    'negative_A': 'A_Host_Random',   'negative_B': 'B_Conserved_Traps', 'negative_C': 'C_EVE_Transposon',
    'host_random': 'A_Host_Random',   'Host_Random': 'A_Host_Random',
    'conserved_trap': 'B_Conserved_Traps', 'Conserved_Traps': 'B_Conserved_Traps',
    'EVE': 'C_EVE_Transposon',       'EVE_Transposon': 'C_EVE_Transposon',
}
# 负样本类型关键词 (用于从 type 列中识别负样本子类型)
NEG_TYPE_KEYWORDS = ['negative', 'host_random', 'Host_Random', 'conserved', 'Conserved', 'EVE', 'EVE_Transposon']

def extract_fp_breakdown(labels_pd, y_pred, tool_label):
    """
    提取假阳性来源构成 — 按序列类型 (type) 分类, 同时计算假阳性率 (FP率 = FP / 该类型总数)
    兼容旧格式 (negative_A/B/C) 和 MASTER 格式 (host_random/conserved_trap/EVE)
    """
    fp_mask = (labels_pd['label'] != 'positive') & (y_pred.astype(bool))
    fp_df = labels_pd[fp_mask]

    result = {'Tool': tool_label, 'Total_FP': int(fp_mask.sum())}
    # 自动检测负样本子类型: 优先 MASTER 关键词, 其次 legacy negative_前缀
    all_types = set(labels_pd['type'].unique())
    neg_types = sorted([t for t in all_types if any(kw in t for kw in NEG_TYPE_KEYWORDS)])
    if not neg_types:
        neg_types = sorted([t for t in all_types if t.startswith('negative')])

    for stype in neg_types:
        total_n = int(labels_pd['type'].isin([stype]).sum())
        fp_n = int(fp_df['type'].isin([stype]).sum())
        label = NEG_TYPE_LABELS.get(stype, stype)
        result[f'FP_{label}'] = fp_n
        result[f'Rate_{label}'] = round(fp_n / max(1, total_n), 4) if total_n > 0 else 0.0
    return result

# ── 绘图模块 ────────────────────────────────────────────────

def set_style(font_scale=1.1):
    sns.set_theme(style='whitegrid', font_scale=font_scale)
    plt.rcParams.update({
        'figure.dpi': 150, 'savefig.dpi': 300,
        'font.family': 'sans-serif',
        'axes.titleweight': 'bold'
    })

def plot_pr_curves(pr_data, outpath):
    """单模式PR曲线 (保留兼容)"""
    set_style(0.9)
    fig, ax = plt.subplots(figsize=(10, 8))
    for label, (prec, rec, ap) in pr_data.items():
        base = label.split(' (')[0]
        ax.plot(rec, prec, linewidth=2.5, color=TOOL_COLORS.get(base, '#333'),
                label=f'{label} (AUPRC={ap:.3f})')
    ax.set_xlabel('Recall (Sensitivity)', fontsize=13)
    ax.set_ylabel('Precision (PPV)', fontsize=13)
    ax.set_title('Precision-Recall Curves', fontsize=16, pad=15)
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=10, frameon=True)
    ax.grid(alpha=0.3, linestyle='--')
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(outpath)}")

def plot_pr_curves_3x1(ds_stores, modes_present, outpath):
    """3×1 PR曲线对比: 每列一个过滤模式"""
    set_style(0.75)
    n_modes = len(modes_present)
    fig, axes = plt.subplots(1, n_modes, figsize=(7.5 * n_modes, 7))
    if n_modes == 1: axes = [axes]

    for ax, mode in zip(axes, modes_present):
        pr_data = ds_stores.get(mode, {})
        if not pr_data: continue
        for label, (prec, rec, ap) in pr_data.items():
            base = label.split(' (')[0]
            ax.plot(rec, prec, linewidth=2.5, color=TOOL_COLORS.get(base, '#333'),
                    label=f'{base} (AUPRC={ap:.3f})')
        ax.set_xlabel('Recall', fontsize=12); ax.set_ylabel('Precision', fontsize=12)
        ax.set_title(FILTER_LABELS[mode], fontsize=15, fontweight='bold', pad=12)
        ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.3, linestyle='--')
        ax.legend(loc='lower left', fontsize=7, frameon=True)

    fig.suptitle('Precision-Recall Curves by Filter Mode', fontsize=16, fontweight='bold', y=1.01)
    fig.subplots_adjust(top=0.9, wspace=0.25)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(outpath)}")

def plot_bar_comparison(overall_df, outpath):
    set_style()
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    metrics = ['Precision', 'Recall', 'F1']
    df = overall_df.set_index('Tool')
    x = np.arange(len(df.index))
    for ax, metric in zip(axes, metrics):
        colors = [TOOL_COLORS.get(t.split(' (')[0], '#333') for t in df.index]
        bars = ax.bar(x, df[metric].astype(float), color=colors, edgecolor='white',
                      alpha=0.85, width=0.7)
        ax.set_title(f'{metric} Score', fontsize=16, fontweight='bold', pad=15)
        ax.set_ylim(0, 1.15)
        for bar, val in zip(bars, df[metric].astype(float)):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.02,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=10,
                    fontweight='bold', color='#333')
        ax.set_xticks(x)
        ax.set_xticklabels(df.index, rotation=45, ha='right', fontsize=12)
        ax.grid(alpha=0.3, axis='y', linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.tight_layout(w_pad=3.0)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(outpath)}")

def plot_coverage_recall(cov_df, outpath):
    """覆盖率召回率: 单张图, ● Raw / ◆ Filtered / ■ Strict, 全局颜色锚定"""
    set_style(0.85)
    df = cov_df.copy()
    def _extract_mode(t):
        for m in FILTER_LABELS.values():
            if m in t: return m
        return 'Raw'
    df['Mode_'] = df['Tool'].apply(_extract_mode)
    df['Base_Tool'] = df['Tool'].apply(lambda x: x.split(' (')[0])

    base_tools = [t for t in _ALL_TOOLS_LIST]
    base_tools = [t for t in base_tools if t in df['Base_Tool'].unique()]
    modes_present = sorted(df['Mode_'].unique(), key=lambda x: FILTER_MODES.index('raw') if x == 'Raw'
                           else (FILTER_MODES.index('filter') if 'Filter' in x and 'Strict' not in x else FILTER_MODES.index('strict')))
    cov_levels = sorted(df['Coverage'].unique(), key=lambda x: int(x.replace('%', '')))
    marker_map = {'Raw': 'o', 'UniProt Filtered': 'D', 'Strict Filtered': 's'}
    mode_ls = {'Raw': '-', 'UniProt Filtered': '--', 'Strict Filtered': '-.'}

    fig, ax = plt.subplots(figsize=(14, 7))
    for mode in modes_present:
        for tool in base_tools:
            sub = df[(df['Base_Tool'] == tool) & (df['Mode_'] == mode)].set_index('Coverage')
            if sub.empty: continue
            ax.plot(cov_levels, [sub.loc[c, 'Recall'] if c in sub.index else np.nan for c in cov_levels],
                    marker=marker_map.get(mode, 'o'), linewidth=2.5, markersize=8,
                    linestyle=mode_ls.get(mode, '-'), color=TOOL_COLORS.get(tool, '#333'),
                    alpha=0.85, markeredgecolor='white', markeredgewidth=1.2)
    ax.set_xlabel('Genome Coverage Level', fontsize=13)
    ax.set_ylabel('Recall', fontsize=13)
    ax.set_title('Recall Sensitivity by Viral Genome Coverage', fontsize=16, pad=15)
    ax.set_ylim(-0.02, 1.08); ax.grid(alpha=0.3, linestyle='--')
    from matplotlib.lines import Line2D
    tool_h = [Line2D([0], [0], color=TOOL_COLORS.get(t, '#333'), linewidth=2.5, label=t) for t in base_tools]
    mode_h = [Line2D([0], [0], marker=m, color='grey', linewidth=2, markersize=8, label=f'{mode}')
              for mode, m in marker_map.items() if mode in modes_present]
    l1 = ax.legend(handles=tool_h, title='Tools', fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left', framealpha=0.9)
    ax.add_artist(l1)
    ax.legend(handles=mode_h, title='Mode (marker)', fontsize=8, bbox_to_anchor=(1.02, 0.5), loc='center left', framealpha=0.9)
    plt.tight_layout(rect=[0, 0, 0.78, 1])
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [coverage] {os.path.basename(outpath)}")

def plot_fp_reduction(overall_df, outpath):
    if 'Mode' not in overall_df.columns or 'FP' not in overall_df.columns: return
    set_style(0.9)
    tools = sorted(overall_df['Tool'].apply(lambda x: x.split(' (')[0]).unique())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'X']
    lines, labels = [], []

    for i, tool in enumerate(tools):
        sub = overall_df[overall_df['Tool'].str.startswith(tool)]
        vals, fp_rates, modes_found = [], [], []
        raw_fp = sub[sub['Mode'] == FILTER_LABELS['raw']]['FP'].values
        raw_fp = int(raw_fp[0]) if len(raw_fp) else 1
        for m in FILTER_MODES:
            row = sub[sub['Mode'] == FILTER_LABELS[m]]
            if len(row):
                fp = int(row['FP'].values[0])
                vals.append(fp); fp_rates.append(fp / max(1, raw_fp) * 100)
                modes_found.append(FILTER_LABELS[m])
        if vals:
            color = PALETTE[i % len(PALETTE)]
            marker = markers[i % len(markers)]
            line, = ax1.plot(modes_found, vals, marker=marker, markersize=10,
                             linewidth=2.5, color=color, alpha=0.85,
                             markeredgecolor='white', markeredgewidth=1.5)
            ax2.plot(modes_found, fp_rates, marker=marker, markersize=10,
                     linewidth=2.5, color=color, alpha=0.85,
                     markeredgecolor='white', markeredgewidth=1.5)
            lines.append(line); labels.append(tool)

    # 左: symlog 防 L 型挤压
    ax1.set_ylabel('Absolute False Positives', fontsize=14, fontweight='bold')
    ax1.set_title('(A) FP Absolute Count Drop', fontsize=16, fontweight='bold', pad=15)
    ax1.set_yscale('symlog', linthresh=5)
    from matplotlib.ticker import ScalarFormatter
    ax1.yaxis.set_major_formatter(ScalarFormatter())

    # 右: 相对比例
    ax2.set_ylabel('Remaining FP (%)', fontsize=14, fontweight='bold')
    ax2.set_title('(B) FP Relative Reduction Rate', fontsize=16, fontweight='bold', pad=15)
    ax2.set_ylim(-5, 105)

    for ax in [ax1, ax2]:
        ax.set_xticklabels(modes_found, fontsize=13)
        ax.grid(alpha=0.3, axis='y', linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.2)
        ax.spines['bottom'].set_linewidth(1.2)

    fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, -0.06),
               ncol=5, fontsize=10, frameon=False, title='Tools', title_fontsize=11)
    fig.suptitle('False Positive Reduction by UniProt Filtering', fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout(rect=[0, 0.04, 1, 0.94])
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [modes] {os.path.basename(outpath)}")

def plot_mode_scatter(overall_df, outpath):
    if 'Mode' not in overall_df.columns: return
    try:
        from adjustText import adjust_text
    except ImportError:
        adjust_text = None
        print("  [info] pip install adjustText 可获得防重叠文本标签")

    set_style(0.9)
    fig, ax = plt.subplots(figsize=(10, 8))

    # 1. F1 ISO-等高线背景
    x_vals = np.linspace(0.01, 1.0, 200)
    y_vals = np.linspace(0.01, 1.0, 200)
    X, Y = np.meshgrid(x_vals, y_vals)
    F1_grid = 2 * (X * Y) / (X + Y + 1e-10)
    CS = ax.contour(X, Y, F1_grid, levels=np.arange(0.1, 1.0, 0.1),
                    colors='gray', alpha=0.25, linestyles='dashed', linewidths=0.8)
    ax.clabel(CS, inline=True, fontsize=8, fmt='F1=%.1f')

    # 2. 双重图例编码: 颜色=工具, 形状=模式
    modes = sorted(overall_df['Mode'].unique())
    markers = {'Raw': 'o', 'UniProt Filtered': 's', 'Strict Filtered': '^'}
    texts = []
    tools_seen = set()
    modes_seen = set()

    for _, r in overall_df.iterrows():
        mode = r['Mode']
        tool = r['Tool'].split(' (')[0]
        color = TOOL_COLORS.get(tool, '#333')
        marker = markers.get(mode, 'o')
        ax.scatter(r['Recall'], r['Precision'], s=120, alpha=0.85,
                   marker=marker, color=color, edgecolor='white', linewidth=1.2, zorder=3)
        tools_seen.add(tool)
        modes_seen.add(mode)
        short = tool[:8]
        texts.append(ax.text(r['Recall'], r['Precision'], short,
                             fontsize=8, fontweight='500', color='#333333'))

    # 3. 智能防重叠
    if adjust_text and texts:
        adjust_text(texts,
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.8, alpha=0.6),
                    expand_points=(1.5, 1.5), force_points=0.5)

    ax.set_xlabel('Recall (Sensitivity)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Precision (PPV)', fontsize=13, fontweight='bold')
    ax.set_title('Precision-Recall Trade-off by Filter Strategy\n(Dashed contours = constant F1)', fontsize=14, fontweight='bold', pad=12)
    ax.set_xlim(0, 1.05); ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.15, linestyle='-')
    # 双图例: 左=工具颜色, 右=模式形状
    from matplotlib.lines import Line2D
    tool_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=TOOL_COLORS.get(t, '#333'),
                            markersize=8, label=t) for t in sorted(tools_seen)]
    mode_handles = [Line2D([0], [0], marker=markers.get(m, 'o'), color='w', markerfacecolor='grey',
                            markersize=8, label=m) for m in sorted(modes_seen)]
    l1 = ax.legend(handles=tool_handles, title='Tool', fontsize=7, title_fontsize=8,
                   loc='lower left', framealpha=0.9)
    ax.add_artist(l1)
    ax.legend(handles=mode_handles, title='Mode (shape)', fontsize=7, title_fontsize=8,
              loc='lower right', framealpha=0.9)
    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight'); plt.close(fig)
    print(f"  [modes] {os.path.basename(outpath)}")

def plot_filter_comparison_3x1(overall_df, outpath):
    if 'Mode' not in overall_df.columns: return
    set_style(0.95)
    df = overall_df.copy()
    df['Base_Tool'] = df['Tool'].apply(lambda x: x.split(' (')[0])
    tools = _ALL_TOOLS_LIST
    tools = [t for t in tools if t in df['Base_Tool'].unique()]
    actual_modes = [m for m in FILTER_LABELS.values() if m in df['Mode'].unique()]
    metrics = ['Recall', 'Precision', 'F1']
    fig, axes = plt.subplots(3, 1, figsize=(15, 16))
    x = np.arange(len(tools)); width = 0.75 / len(actual_modes)
    mode_colors = ['#76C8A6', '#8FAADC', '#FDE05A']
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        for i, mode in enumerate(actual_modes):
            md = df[df['Mode'] == mode].set_index('Base_Tool')
            vals = [float(md.loc[t, metric]) if t in md.index else 0 for t in tools]
            ax.bar(x + i * width, vals, width, alpha=0.9, label=mode,
                   color=mode_colors[i % len(mode_colors)], edgecolor='white', linewidth=0.8)
        ax.set_xticks(x + width * (len(actual_modes) - 1) / 2)
        if idx == 2:
            ax.set_xticklabels(tools, rotation=45, ha='right', fontsize=13)
        else:
            ax.set_xticklabels([])
        ax.set_ylabel(metric, fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.2); ax.grid(alpha=0.3, axis='y', linestyle='--')
    # 图例放最下方子图的下面
    axes[-1].legend(fontsize=12, loc='upper center', bbox_to_anchor=(0.5, -0.12),
                    ncol=3, frameon=True)
    fig.suptitle('Metrics Distribution Across Filter Modes', fontsize=17, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0.02, 1, 0.95], h_pad=3.0)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [modes] {os.path.basename(outpath)}")

def plot_radar_charts(overall_df, outpath):
    """3×1 雷达图组合: 每行一个过滤模式, P/R/F1 三轴"""
    if 'Mode' not in overall_df.columns: return
    set_style(0.8)
    df = overall_df.copy()
    df['Base_Tool'] = df['Tool'].apply(lambda x: x.split(' (')[0])
    actual_modes = [m for m in FILTER_LABELS.values() if m in df['Mode'].unique()]
    if len(actual_modes) < 1: return

    md_first = df[df['Mode'] == actual_modes[0]].set_index('Base_Tool')
    tools = sorted(md_first.index.tolist())
    if 'Ensemble' in tools: tools.remove('Ensemble'); tools.append('Ensemble')

    num_vars = len(tools)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]
    metrics = ['Precision', 'Recall', 'F1']
    metric_colors = {'Precision': '#4DBBD5', 'Recall': '#00A087', 'F1': '#E64B35'}

    n_modes = len(actual_modes)
    fig, axes = plt.subplots(1, n_modes, figsize=(8 * n_modes, 8),
                              subplot_kw=dict(polar=True))
    if n_modes == 1:
        axes = [axes]

    for ax, mode in zip(axes, actual_modes):
        md = df[df['Mode'] == mode].set_index('Base_Tool')
        for metric in metrics:
            values = [float(md.loc[t, metric]) if t in md.index else 0 for t in tools]
            values += values[:1]
            color = metric_colors[metric]
            ax.plot(angles, values, color=color, linewidth=2, label=metric)
            ax.fill(angles, values, color=color, alpha=0.1)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(tools, fontsize=10)
        ax.tick_params(axis='x', pad=18)
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], color='grey', size=8)
        ax.set_title('')  # 清除默认标题
        ax.text(0.5, -0.15, mode, transform=ax.transAxes, ha='center', va='top',
                fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', bbox_to_anchor=(1.32, 1.08), fontsize=10, framealpha=0.9)

    fig.suptitle('Performance Radar — Precision / Recall / F1', fontsize=16, fontweight='bold', y=0.98)
    fig.subplots_adjust(top=0.88, wspace=0.3, bottom=0.12)
    plt.tight_layout(pad=3.0)
    out_file = os.path.join(outpath, 'Fig_Radar_3x1.png')
    fig.savefig(out_file, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(out_file)}")

def _to_upset_contents(y_preds_dict):
    """去后缀, 构建 upsetplot from_contents 需要的 dict (np.int64 → Python int)"""
    contents = {}
    for t, preds in y_preds_dict.items():
        base = t.split(' (')[0] if ' (' in t else t
        if base == 'Ensemble': continue
        contents[base] = set(int(x) for x in np.where(preds)[0])
    return contents

def _build_upset_from_preds(y_preds_dict):
    """从预测dict构建UpSet数据: tools, set_sizes, combo_tools, combo_counts, matrix"""
    tools = []
    seen = set()
    for t in y_preds_dict.keys():
        base = t.split(' (')[0] if ' (' in t else t
        if base != 'Ensemble' and base not in seen:
            seen.add(base); tools.append(base)
    if not tools: return None
    preds = {t: y_preds_dict[t].astype(bool) for t in tools if t in y_preds_dict}
    # 补充: 从带后缀的key提取
    for t in list(tools):
        if t not in preds:
            for k, v in y_preds_dict.items():
                if k.split(' (')[0] == t:
                    preds[t] = v.astype(bool); break
    df = pd.DataFrame(preds)
    mask = df.any(axis=1); df_sub = df[mask]
    groups = df_sub.groupby(tools).size().sort_values(ascending=False).head(25)
    combo_tools, combo_counts = [], []
    for idx, count in groups.items():
        active = [t for t, v in zip(tools, idx) if v]
        if active: combo_tools.append(active); combo_counts.append(count)
    if not combo_tools: return None
    set_sizes = [int(np.sum(preds[t])) for t in tools]
    n = len(combo_tools)
    matrix = np.zeros((len(tools), n))
    for j, ct in enumerate(combo_tools):
        for t in ct: matrix[tools.index(t), j] = 1
    return tools, set_sizes, combo_tools, combo_counts, matrix

def _draw_upset_cell(fig, gs_cell, tools, set_sizes, combo_tools, combo_counts, matrix, title, colors_dict):
    """在GridSpec格子中绘制单个UpSet (左set + 右上柱 + 右下点阵)"""
    gs = gs_cell.subgridspec(2, 2, height_ratios=[1, 2.5], width_ratios=[1, 4], hspace=0.05, wspace=0.05)
    ax_set = fig.add_subplot(gs[:, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)

    # 左: set size 水平柱
    colors = [colors_dict.get(t, '#333') for t in tools]
    ax_set.barh(range(len(tools)), set_sizes, color=colors, edgecolor='white', height=0.7)
    ax_set.set_xlim(0, max(set_sizes)*1.2); ax_set.invert_yaxis()
    ax_set.set_yticks(range(len(tools)))
    ax_set.set_yticklabels(tools, fontsize=8); ax_set.tick_params(axis='y', length=0)
    ax_set.spines['top'].set_visible(False); ax_set.spines['right'].set_visible(False)
    for i, v in enumerate(set_sizes):
        ax_set.text(v+max(set_sizes)*0.02, i, str(v), va='center', fontsize=7)

    # 右上: 交集柱状图
    n = len(combo_counts)
    ccolors = []
    for xi in range(n):
        a = combo_tools[xi]
        ccolors.append(colors_dict.get(a[0], '#3C5488') if len(a)==1 else '#3C5488')
    x = np.arange(n)
    ax_bar.bar(x, combo_counts, color=ccolors, edgecolor='white', width=0.7)
    for i, v in enumerate(combo_counts):
        ax_bar.text(i, v+max(combo_counts)*0.03, str(v), ha='center', fontsize=7)
    ax_bar.grid(alpha=0.3, axis='y', linestyle='--')
    ax_bar.spines['top'].set_visible(False); ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['bottom'].set_visible(False); ax_bar.tick_params(labelbottom=False)

    # 右下: 点阵矩阵
    for yi in range(len(tools)):
        ax_mat.axhline(yi, color='gray', alpha=0.10, linestyle='-', zorder=1)
    for xi in range(n):
        ay = [yi for yi in range(len(tools)) if matrix[yi, xi]==1]
        if ay:
            c = ccolors[xi]
            ax_mat.plot([xi, xi], [min(ay), max(ay)], color=c, linewidth=2.5, zorder=2)
            ax_mat.scatter([xi]*len(ay), ay, color=c, s=100, zorder=3)
        iy = [yi for yi in range(len(tools)) if matrix[yi, xi]==0]
        if iy: ax_mat.scatter([xi]*len(iy), iy, color='#E8E8E8', s=30, zorder=2)
    ax_mat.set_xticks(range(n))
    ax_mat.set_xticklabels([f'C{i+1}' for i in range(n)], rotation=45, ha='right', fontsize=7)
    ax_mat.spines['top'].set_visible(False); ax_mat.spines['right'].set_visible(False)
    ax_mat.spines['left'].set_visible(False); ax_mat.tick_params(axis='y', left=False, labelleft=False)
    ax_set.set_title(title, fontsize=11, fontweight='bold')

def plot_overlap_3x1(y_preds_per_mode, modes_present, outpath):
    """1×3 UpSet: matplotlib手动绘制, TOOL_COLORS着色, 3列并排"""
    set_style(0.8)
    n_modes = len(modes_present)
    fig = plt.figure(figsize=(8 * n_modes, 7))
    gs = fig.add_gridspec(1, n_modes, wspace=0.08, left=0.02, right=0.98, top=0.90, bottom=0.08)

    for idx, mode in enumerate(modes_present):
        yp = y_preds_per_mode.get(mode, {})
        if not yp: continue
        data = _build_upset_from_preds(yp)
        if data is None: continue
        _draw_upset_cell(fig, gs[idx], *data, FILTER_LABELS[mode], TOOL_COLORS)

    fig.suptitle('Detection Overlap — UpSet by Filter Mode', fontsize=15, fontweight='bold', y=0.97)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overlap] {os.path.basename(outpath)}")

def plot_overlap(y_preds_dict, outpath):
    """单模式UpSet — upsetplot原生 + 色彩"""
    contents = _to_upset_contents(y_preds_dict)
    total_elements = sum(len(v) for v in contents.values())
    if len(contents) < 2 or total_elements == 0: return
    try:
        from upsetplot import from_contents, UpSet
        data = from_contents(contents)
        fig = plt.figure(figsize=(14, 8))
        upset = UpSet(data, subset_size='count',
                      facecolor='#3C5488',
                      other_dots_color='#CCCCCC',
                      shading_color='#F0F2F6',
                      min_subset_size=max(1, int(total_elements * 0.001)),
                      show_counts='%d', sort_by='cardinality')
        upset.plot(fig=fig)
        fig.suptitle('Detection Overlap — UpSet Plot', fontsize=16, fontweight='bold', y=1.02)
        fig.savefig(outpath, bbox_inches='tight', dpi=300); plt.close(fig)
        print(f"  [overlap] {os.path.basename(outpath)}")
    except Exception as e:
        print(f"  [overlap] upsetplot error: {e}")

def _plot_overlap_simple(y_preds_dict, outpath):
    """降级方案: 共识度分布柱状图"""
    set_style(0.9)
    tools_bar = [t for t in y_preds_dict.keys() if t != 'Ensemble']
    if len(tools_bar) < 2: return
    votes = np.sum(np.column_stack([y_preds_dict[t] for t in tools_bar]), axis=1)
    counts = [int(np.sum(votes == i)) for i in range(1, len(tools_bar) + 1)]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(range(1, len(tools_bar) + 1), counts,
           color=plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(tools_bar))), edgecolor='white')
    ax.set_xlabel('Tools detecting consensus', fontsize=12)
    ax.set_ylabel('Sequences', fontsize=12)
    ax.set_title('Detection Consensus', fontsize=15, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    for i, v in enumerate(counts):
        if v > 0: ax.text(i + 1, v + max(counts) * 0.02, str(v), ha='center', fontsize=11)
    plt.tight_layout()
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overlap] {os.path.basename(outpath)}")

def plot_combined_figure(ds_store, overall_basic, cov_df, outpath):
    """ABCD 四合一组合图: 超大画布 + 图例外置 + 斑马线表格"""
    set_style(0.8)
    fig = plt.figure(figsize=(22, 16))

    # (A) PR Curves
    ax1 = fig.add_subplot(2, 2, 1)
    for label, (prec, rec, ap) in ds_store.items():
        base = label.split(' (')[0]
        ax1.plot(rec, prec, linewidth=2.5, color=TOOL_COLORS.get(base, '#333'),
                 label=f'{label} (AUPRC={ap:.3f})')
    ax1.set_xlabel('Recall', fontsize=12); ax1.set_ylabel('Precision', fontsize=12)
    ax1.set_title('(A) PR Curves', fontsize=16, fontweight='bold')
    ax1.set_xlim(0, 1.05); ax1.set_ylim(0, 1.05); ax1.grid(alpha=0.3, linestyle='--')
    ax1.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=True)

    # (B) Method Comparison
    ax2 = fig.add_subplot(2, 2, 2)
    df_bar = overall_basic.set_index('Tool')
    x = np.arange(len(df_bar)); w = 0.25
    metric_colors_b = ['#8FAADC', '#76C8A6', '#FDE05A']
    for idx, metric in enumerate(['Precision', 'Recall', 'F1']):
        ax2.bar(x + idx * w, df_bar[metric].astype(float), w, alpha=0.9,
                color=metric_colors_b[idx], label=metric)
    ax2.set_xticks(x + w)
    ax2.set_xticklabels(df_bar.index, rotation=35, ha='right', fontsize=11)
    ax2.set_ylim(0, 1.25); ax2.set_ylabel('Score', fontsize=12)
    ax2.set_title('(B) Baseline Metrics', fontsize=16, fontweight='bold')
    ax2.legend(loc='upper right', ncol=3, fontsize=10)
    ax2.grid(alpha=0.3, axis='y', linestyle='--')

    # (C) Recall by Coverage
    ax3 = fig.add_subplot(2, 2, 3)
    if len(cov_df):
        cov_levels = sorted(cov_df['Coverage'].unique(), key=lambda x: int(x.replace('%', '')))
        for tool in cov_df['Tool'].unique():
            sub = cov_df[cov_df['Tool'] == tool].set_index('Coverage')
            base = tool.split(' (')[0]
            ax3.plot(cov_levels, [sub.loc[c, 'Recall'] for c in cov_levels],
                     marker='o', linewidth=2.5, markersize=7,
                     color=TOOL_COLORS.get(base, '#333'), label=tool)
        ax3.set_xlabel('Coverage Level', fontsize=12); ax3.set_ylabel('Recall', fontsize=12)
        ax3.set_title('(C) Recall Sensitivity by Coverage', fontsize=16, fontweight='bold')
        ax3.set_ylim(0, 1.1); ax3.grid(alpha=0.3, linestyle='--')
        ax3.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=True)

    # (D) Data Table
    ax4 = fig.add_subplot(2, 2, 4); ax4.axis('off')
    tbl_cols = [c for c in ['Tool', 'Precision', 'Recall', 'F1', 'AUPRC'] if c in overall_basic.columns]
    tbl_data = overall_basic[tbl_cols].round(4)
    tbl = ax4.table(cellText=tbl_data.values, colLabels=tbl_data.columns, cellLoc='center', loc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(11); tbl.scale(1.1, 1.8)
    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#3C5488')
        elif row % 2 == 0:
            cell.set_facecolor('#F3F6F9')  # 斑马纹
    ax4.set_title('(D) Quantitative Summary', fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout(w_pad=8.0, h_pad=4.0)
    fig.savefig(outpath, bbox_inches='tight'); plt.close(fig)
    print(f"  [overall] {os.path.basename(outpath)}")

def plot_fp_breakdown_dual(fp_df, outpath, mode_label='raw'):
    """单模式假阳性双拼图: 左=FP率(归一化) + 右=绝对数量(symlog)"""
    if fp_df is None or (hasattr(fp_df, 'empty') and fp_df.empty) or 'Tool' not in fp_df.columns:
        return
    set_style()
    df_raw = fp_df.copy().set_index('Tool').fillna(0)
    fp_cols = sorted([c for c in df_raw.columns if c.startswith('FP_')])
    rate_cols = sorted([c for c in df_raw.columns if c.startswith('Rate_')])
    if not fp_cols: return

    df_fp = df_raw[fp_cols]
    df_fp.columns = [c.replace('FP_', '') for c in fp_cols]
    df_rate = df_raw[rate_cols] if rate_cols else pd.DataFrame()
    if not df_rate.empty:
        df_rate.columns = [c.replace('Rate_', '') for c in rate_cols]

    colors = PALETTE[:len(df_fp.columns)]

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # 左: 假阳性率 (FP / 该类型总数, 消除样本量差异)
    if not df_rate.empty:
        df_rate.plot(kind='bar', stacked=False, color=colors, ax=axes[0],
                     edgecolor='white', linewidth=0.5, width=0.7)
    else:
        df_fp.plot(kind='bar', stacked=True, color=colors, ax=axes[0],
                   edgecolor='white', linewidth=0.5, width=0.75)
    axes[0].set_title('(A) FP Rate (FP / Category Size)', fontsize=14, fontweight='bold')
    axes[0].set_ylabel('False Positive Rate', fontsize=12)
    axes[0].set_xticklabels(df_fp.index, rotation=35, ha='right')
    axes[0].legend(title='FP Source', bbox_to_anchor=(1.05, 1), loc='upper left',
                   fontsize=9, framealpha=0.9)
    axes[0].grid(alpha=0.3, axis='y', linestyle='--')

    # 右: 绝对数量 (symlog)
    df_fp.plot(kind='bar', stacked=True, color=colors, ax=axes[1],
               edgecolor='white', linewidth=0.5, width=0.75)
    axes[1].set_title('(B) FP Absolute Count (SymLog)', fontsize=14, fontweight='bold')
    axes[1].set_ylabel('Number of False Positives', fontsize=12)
    axes[1].set_yscale('symlog', linthresh=10)
    axes[1].set_xticklabels(df_fp.index, rotation=35, ha='right')
    axes[1].get_legend().remove()
    axes[1].grid(alpha=0.3, axis='y', linestyle='--')

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    out_file = os.path.join(outpath, f'Fig_FP_Breakdown_Dual_{mode_label.replace(" ", "_")}.png')
    fig.savefig(out_file, bbox_inches='tight'); plt.close(fig)
    print(f"  [entity] {os.path.basename(out_file)}")


def plot_fp_breakdown_dual_3x1(all_fp_dfs, modes_present, outpath):
    """3×1 假阳性双拼图: 每行一个过滤模式, 左=FP率(归一化) + 右=绝对数量(symlog)"""
    if not all_fp_dfs or not modes_present: return
    set_style(0.8)
    n_modes = len(modes_present)

    fp_cols_all = set()
    rate_cols_all = set()
    for fp_df in all_fp_dfs:
        if fp_df is None or (hasattr(fp_df, 'empty') and fp_df.empty) or 'Tool' not in fp_df.columns:
            continue
        fp_cols_all.update([c for c in fp_df.columns if c.startswith('FP_')])
        rate_cols_all.update([c for c in fp_df.columns if c.startswith('Rate_')])
    fp_cols_all = sorted(fp_cols_all)
    rate_cols_all = sorted(rate_cols_all)
    if not fp_cols_all: return
    colors = PALETTE[:len(fp_cols_all)]

    fig, axes = plt.subplots(n_modes, 2, figsize=(18, 6 * n_modes))
    if n_modes == 1: axes = [axes]

    for row_idx, (mode, fp_df) in enumerate(zip(modes_present, all_fp_dfs)):
        if fp_df is None or (hasattr(fp_df, 'empty') and fp_df.empty):
            axes[row_idx][0].set_visible(False); axes[row_idx][1].set_visible(False); continue

        df_raw = fp_df.copy().set_index('Tool').fillna(0)
        for c in fp_cols_all:
            if c not in df_raw.columns: df_raw[c] = 0
        for c in rate_cols_all:
            if c not in df_raw.columns: df_raw[c] = 0.0

        df_fp = df_raw[fp_cols_all]
        df_fp.columns = [c.replace('FP_', '') for c in fp_cols_all]
        df_rate = df_raw[rate_cols_all] if rate_cols_all else pd.DataFrame()
        if not df_rate.empty:
            df_rate.columns = [c.replace('Rate_', '') for c in rate_cols_all]

        ax1, ax2 = axes[row_idx]
        if not df_rate.empty:
            df_rate.plot(kind='bar', stacked=False, color=colors, ax=ax1, edgecolor='white', linewidth=0.5, width=0.7)
        else:
            df_fp.plot(kind='bar', stacked=True, color=colors, ax=ax1, edgecolor='white', linewidth=0.5, width=0.75)
        ax1.set_ylabel(f'{FILTER_LABELS[mode]}', fontsize=11, fontweight='bold')
        ax1.set_xticklabels(df_fp.index, rotation=30, ha='right', fontsize=9)
        ax1.grid(alpha=0.3, axis='y', linestyle='--')
        df_fp.plot(kind='bar', stacked=True, color=colors, ax=ax2, edgecolor='white', linewidth=0.5, width=0.75)
        ax2.set_ylabel(f'{FILTER_LABELS[mode]}', fontsize=11, fontweight='bold')
        ax2.set_yscale('symlog', linthresh=10)
        ax2.set_xticklabels(df_fp.index, rotation=30, ha='right', fontsize=9)
        ax2.grid(alpha=0.3, axis='y', linestyle='--')

    # 全局图例 (来自第一行)
    from matplotlib.patches import Patch
    fp_labels_clean = [c.replace('FP_', '') for c in fp_cols_all]
    legend_patches = [Patch(color=c, label=l) for c, l in zip(colors, fp_labels_clean)]
    fig.legend(handles=legend_patches, title='FP Source', loc='upper right',
               bbox_to_anchor=(0.99, 0.99), fontsize=9, title_fontsize=10, framealpha=0.9)

    fig.suptitle('False Positive Breakdown — Rate (left) vs Count (right)', fontsize=15, fontweight='bold', y=1.01)
    fig.subplots_adjust(top=0.94, hspace=0.55, wspace=0.25, right=0.88)
    out_file = os.path.join(outpath, 'Fig_FP_Breakdown_Dual_3x1.png')
    fig.savefig(out_file, bbox_inches='tight'); plt.close(fig)
    print(f"  [entity] {os.path.basename(out_file)}")

def plot_entity_dumbbell_3x1(all_entity_dfs, modes_present, outpath):
    """3×1 哑铃图合并: 每行一个过滤模式"""
    if not all_entity_dfs or not modes_present: return
    set_style(0.85)
    n_modes = len(modes_present)
    fig, axes = plt.subplots(1, n_modes, figsize=(6 * n_modes, 7))
    if n_modes == 1: axes = [axes]

    # 取第一个非空 DataFrame 作为共享 Y 轴的工具列表 (提取基工具名)
    shared_tools = None
    for ent_df in all_entity_dfs:
        if ent_df is not None and not (hasattr(ent_df, 'empty') and ent_df.empty):
            ent_df = ent_df.copy()
            ent_df['Base_Tool'] = ent_df['Tool'].apply(lambda x: x.split(' (')[0])
            shared_tools = ent_df.sort_values('Lax_Entity_Recall', ascending=True)['Base_Tool'].tolist()
            break
    if not shared_tools: return

    # 只用左轴显示工具名，其余不显示Y标签
    for ax_idx, (ax, mode, ent_df) in enumerate(zip(axes, modes_present, all_entity_dfs)):
        if ent_df is None or (hasattr(ent_df, 'empty') and ent_df.empty):
            ax.set_visible(False); continue
        df = ent_df.copy()
        df['Base_Tool'] = df['Tool'].apply(lambda x: x.split(' (')[0])
        present_tools = [t for t in shared_tools if t in df['Base_Tool'].values]
        if not present_tools: continue
        df = df.set_index('Base_Tool').loc[present_tools].reset_index()
        y_range = range(1, len(shared_tools) + 1)
        ax.hlines(y=y_range, xmin=df['Strict_Entity_Recall'], xmax=df['Lax_Entity_Recall'],
                  color='grey', alpha=0.5, linewidth=3, zorder=1)
        ax.scatter(df['Strict_Entity_Recall'], y_range, color='#E64B35', s=100,
                   label='Strict', zorder=2, edgecolors='white', linewidth=0.5)
        ax.scatter(df['Lax_Entity_Recall'], y_range, color='#4DBBD5', s=100,
                   label='Lax', zorder=3, edgecolors='white', linewidth=0.5)
        for y, strict, lax in zip(y_range, df['Strict_Entity_Recall'], df['Lax_Entity_Recall']):
            ax.text(strict - 0.03, y, f'{strict:.2f}', va='center', ha='right',
                    fontsize=8, color='#E64B35', fontweight='bold')
            ax.text(lax + 0.03, y, f'{lax:.2f}', va='center', ha='left',
                    fontsize=8, color='#4DBBD5', fontweight='bold')
        # 仅左轴显示工具名
        if ax_idx == 0:
            ax.set_yticks(y_range)
            ax.set_yticklabels(shared_tools, fontsize=10)
        else:
            ax.set_yticks([])
        # 模式名放图下方
        ax.set_xlabel(FILTER_LABELS[mode], fontsize=13, fontweight='bold')
        ax.set_xlim(-0.05, 1.1)
        ax.grid(alpha=0.3, axis='x', linestyle='--')
        if ax_idx == n_modes - 1:
            ax.legend(loc='lower right', fontsize=9, framealpha=0.9)

    fig.suptitle('Segmented Virus Recovery Gap', fontsize=15, fontweight='bold', y=0.99)
    fig.subplots_adjust(top=0.92, wspace=0.08)
    out_file = os.path.join(outpath, 'Fig_Entity_Dumbbell_3x1.png')
    fig.savefig(out_file, bbox_inches='tight'); plt.close(fig)
    print(f"  [entity] {os.path.basename(out_file)}")

def plot_entity_dumbbell(entity_df, outpath, mode_label):
    """单模式哑铃图 (保留兼容)"""
    plot_entity_dumbbell_3x1([entity_df], [mode_label], outpath)

def plot_combinations(combo_df, outpath):
    """UpSet风格矩阵图: 上方F1柱状图 + 下方工具参与点阵"""
    if combo_df.empty: return
    set_style(0.85)
    top20 = combo_df.sort_values('F1', ascending=False).head(20).reset_index(drop=True)
    f1_scores = top20['F1'].astype(float).values
    strategies = top20['Strategy'].values

    # 确定所有出现过的工具
    all_tools = _ALL_TOOLS_LIST  # 来自全局配置
    all_tools = [t for t in all_tools if t != 'Ensemble']

    # 构建 0/1 矩阵: 行=工具, 列=组合
    matrix = np.zeros((len(all_tools), len(top20)))
    for i, row in top20.iterrows():
        subset_str = row['Subset']
        for j, tool in enumerate(all_tools):
            if tool.lower() in subset_str.lower():
                matrix[j, i] = 1

    # 上下布局: 30% 柱状图 + 70% 点阵矩阵
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 2.5], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)

    # ── 上方: F1柱状图 ──
    strat_colors = {'Union': '#E64B35', 'Majority': '#4DBBD5', 'Intersection': '#00A087'}
    colors = [strat_colors.get(s, '#333') for s in strategies]
    x = np.arange(len(top20))
    bars = ax_bar.bar(x, f1_scores, color=colors, edgecolor='white', width=0.6)
    min_f1 = f1_scores.min()
    ax_bar.set_ylim(max(0, min_f1 - 0.02), f1_scores.max() + 0.005)
    ax_bar.set_ylabel('F1 Score', fontsize=13, fontweight='bold')
    for bar, val in zip(bars, f1_scores):
        ax_bar.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                    f'{val:.4f}', ha='center', va='bottom', fontsize=10, rotation=45)
    ax_bar.grid(alpha=0.3, axis='y', linestyle='--')
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['bottom'].set_visible(False)
    ax_bar.tick_params(labelbottom=False, bottom=False)

    # ── 下方: UpSet点阵矩阵 ──
    y = np.arange(len(all_tools))
    for yi in y:
        ax_matrix.axhline(yi, color='gray', alpha=0.12, linestyle='-', zorder=1)
    for xi in x:
        active_y = [yi for yi, active in enumerate(matrix[:, xi]) if active == 1]
        if active_y:
            ax_matrix.plot([xi, xi], [min(active_y), max(active_y)], color='#555555', linewidth=2.5, zorder=2)
            ax_matrix.scatter([xi] * len(active_y), active_y, color='#3C5488', s=150, zorder=3)
        inactive_y = [yi for yi, active in enumerate(matrix[:, xi]) if active == 0]
        if inactive_y:
            ax_matrix.scatter([xi] * len(inactive_y), inactive_y, color='#DDDDDD', s=50, zorder=2)

    ax_matrix.set_yticks(y)
    ax_matrix.set_yticklabels(all_tools, fontsize=12)
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels([f'R{i+1}' for i in x], rotation=45, ha='right', fontsize=11)
    ax_matrix.invert_yaxis()
    ax_matrix.spines['top'].set_visible(False)
    ax_matrix.spines['right'].set_visible(False)
    ax_matrix.spines['bottom'].set_visible(False)
    ax_matrix.spines['left'].set_visible(False)
    ax_matrix.tick_params(axis='y', length=0)

    # ── 标题在上, 图例紧贴标题下方 ──
    import matplotlib.patches as mpatches
    legend_patches = [mpatches.Patch(color=c, label=s) for s, c in strat_colors.items()]
    fig.suptitle('Top 20 Ensemble Combinations — Composition & Performance', fontsize=16, fontweight='bold', y=0.99)
    fig.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, 0.94), ncol=3,
               fontsize=10, frameon=False, columnspacing=1.5)

    fig.subplots_adjust(top=0.90, bottom=0.08, hspace=0.05)
    fig.savefig(os.path.join(outpath, 'Fig_Combinations_UpSet.png'), dpi=300, bbox_inches='tight')

    plt.close(fig); print(f"  [combinations] Fig_Combinations_UpSet.png")

def plot_marginal_benefit(combo_df, outpath):
    """边际收益散点图: 工具数量 vs F1, 展示投入产出比"""
    if combo_df.empty or 'N_Tools' not in combo_df.columns: return
    set_style()
    fig, ax = plt.subplots(figsize=(12, 7))
    strat_styles = {
        'Union': {'color': '#E64B35', 'marker': 'o'},
        'Majority': {'color': '#4DBBD5', 'marker': 's'},
        'Intersection': {'color': '#00A087', 'marker': '^'}
    }
    np.random.seed(42)
    for strat, style in strat_styles.items():
        subset = combo_df[combo_df['Strategy'] == strat]
        if subset.empty: continue
        jitter = np.random.uniform(-0.15, 0.15, size=len(subset))
        ax.scatter(subset['N_Tools'] + jitter, subset['F1'].astype(float),
                   color=style['color'], marker=style['marker'],
                   s=60, alpha=0.55, edgecolor='white', linewidth=0.5, label=strat)

    combo_df['F1'] = combo_df['F1'].astype(float)
    max_f1 = combo_df.groupby('N_Tools')['F1'].max().sort_index()
    ax.plot(max_f1.index, max_f1.values, color='#333333', linewidth=2.5, linestyle='--', zorder=4,
            label='Maximum F1 Envelope (Diminishing Returns)')
    for n in max_f1.index:
        if n in [1, 2, 3, 9] or max_f1[n] == max_f1.max():
            ax.scatter(n, max_f1[n], color='#333333', s=100, zorder=5)
            ax.text(n, max_f1[n] + 0.01, f'{max_f1[n]:.3f}', ha='center', va='bottom', fontweight='bold', fontsize=9)

    ax.set_xticks(range(1, 10))
    ax.set_xlabel('Integration Complexity (Number of Tools)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Performance (F1 Score)', fontsize=13, fontweight='bold')
    ax.set_title('Marginal Benefit Analysis of Tool Ensembles', fontsize=16, fontweight='bold', pad=15)
    ax.grid(alpha=0.3, linestyle='--')
    ax.legend(title='', loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=3, fontsize=10, frameon=False)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    min_f1_all = combo_df['F1'].astype(float).min()
    ax.set_ylim(max(0, min_f1_all - 0.05), max_f1.max() + 0.05)
    fig.subplots_adjust(bottom=0.12)
    fig.savefig(os.path.join(outpath, 'Fig_Combinations_Marginal_Benefit.png'), dpi=300, bbox_inches='tight')
    plt.close(fig); print(f"  [combinations] Fig_Combinations_Marginal_Benefit.png")

def plot_synergy_heatmap(combo_df, outpath):
    """组合协同热力图: 左侧工具包含矩阵 + 右侧F1柱状图"""
    if combo_df.empty: return
    set_style(0.85)
    top20 = combo_df.sort_values('F1', ascending=False).head(20).reset_index(drop=True)
    f1_scores = top20['F1'].astype(float).values
    all_tools = [t for t in _ALL_TOOLS_LIST if t != 'Ensemble']

    matrix = np.zeros((len(top20), len(all_tools)))
    for i, row in top20.iterrows():
        for j, tool in enumerate(all_tools):
            if tool.lower() in row['Subset'].lower():
                matrix[i, j] = 1

    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1], wspace=0.05)
    ax_heat = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1], sharey=ax_heat)

    cmap = sns.color_palette(["#F0F2F6", "#3C5488"])
    sns.heatmap(matrix, cmap=cmap, cbar=False, ax=ax_heat, linewidths=1.5, linecolor='white')
    ax_heat.set_xticks(np.arange(len(all_tools)) + 0.5)
    ax_heat.set_xticklabels(all_tools, rotation=45, ha='right', fontsize=11, fontweight='bold')
    y_labels = [f"R{i+1} [{row['Strategy'][:3]}]" for i, row in top20.iterrows()]
    ax_heat.set_yticks(np.arange(len(top20)) + 0.5)
    ax_heat.set_yticklabels(y_labels, rotation=0, fontsize=9)
    ax_heat.set_title('Tool Presence Matrix', fontsize=14, fontweight='bold', pad=10)

    y_pos = np.arange(len(top20)) + 0.5
    norm = plt.Normalize(f1_scores.min() - 0.01, f1_scores.max())
    colors = plt.cm.YlOrRd(norm(f1_scores))
    ax_bar.barh(y_pos, f1_scores, color=colors, height=0.7, edgecolor='white')
    min_f1 = f1_scores.min()
    ax_bar.set_xlim(max(0, min_f1 - 0.02), f1_scores.max() + 0.005)
    for y, val in zip(y_pos, f1_scores):
        ax_bar.text(val + 0.002, y, f'{val:.4f}', va='center', fontsize=9, fontweight='bold')
    ax_bar.set_title('F1 Score', fontsize=14, fontweight='bold', pad=10)
    ax_bar.spines['top'].set_visible(False); ax_bar.spines['right'].set_visible(False); ax_bar.spines['left'].set_visible(False)
    ax_bar.tick_params(axis='y', left=False, labelleft=False)
    ax_bar.grid(alpha=0.3, axis='x', linestyle='--')

    fig.suptitle('Synergy Map of Optimal Tool Combinations', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(os.path.join(outpath, 'Fig_Combinations_Synergy_Heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close(fig); print(f"  [combinations] Fig_Combinations_Synergy_Heatmap.png")


# ── 组合分析 ────────────────────────────────────────────────

def run_combination_analysis(y_preds_raw, labels_pd, y_true_all, outpath, tool_names_sorted):
    """穷举组合 (2^n-1 种) × 3 投票策略 + P1-P5 预设分组"""
    if not y_preds_raw:
        return pd.DataFrame()
    y_pred_matrix = np.column_stack([y_preds_raw[t] for t in tool_names_sorted])
    n_tools = len(tool_names_sorted)
    combo_results = []
    for r in range(1, n_tools + 1):
        for subset in combinations(range(n_tools), r):
            names = '+'.join([tool_names_sorted[i] for i in subset])
            n_sel = len(subset)
            votes = np.sum(y_pred_matrix[:, subset], axis=1)
            for s_label, pred in [('Union', votes >= 1),
                                   ('Majority', votes >= max(1, n_sel // 2 + 1)),
                                   ('Intersection', votes >= n_sel)]:
                m = compute_metrics(y_true_all, pred.astype(int), f'{names} ({s_label})')
                ap, _, _, _ = compute_auprc(y_true_all, pred.astype(float))
                m['Strategy'] = s_label; m['Subset'] = names; m['N_Tools'] = n_sel
                m['AUPRC'] = round(ap, 4); combo_results.append(m)
    combo_df = pd.DataFrame(combo_results)
    if not combo_df.empty:
        combo_df.to_csv(os.path.join(outpath, 'identification_combinations.tsv'), sep='\t', index=False)
    # P1-P5 预设
    pre_rows = []
    for combo_name, tool_subset in PREDEFINED_COMBOS.items():
        indices = [i for i, t in enumerate(tool_names_sorted) if any(t.startswith(b) for b in tool_subset)]
        if not indices: continue
        votes = np.sum(y_pred_matrix[:, indices], axis=1)
        n_sel = len(indices)
        for s_label, pred in [('Union', votes >= 1),
                               ('Majority', votes >= max(1, n_sel // 2 + 1)),
                               ('Intersection', votes >= n_sel)]:
            m = compute_metrics(y_true_all, pred.astype(int), f'{combo_name} ({s_label})')
            ap, _, _, _ = compute_auprc(y_true_all, pred.astype(float))
            m['Strategy'] = s_label; m['Group'] = combo_name; m['N_Tools'] = n_sel
            m['AUPRC'] = round(ap, 4); pre_rows.append(m)
    pre_df = pd.DataFrame(pre_rows) if pre_rows else pd.DataFrame()
    if not pre_df.empty:
        pre_df.to_csv(os.path.join(outpath, 'identification_predefined_combos.tsv'), sep='\t', index=False)
        print(f"\n  P1-P5 预设组合报告:")
        cols = [c for c in ['Group', 'Strategy', 'Precision', 'Recall', 'F1', 'AUPRC'] if c in pre_df.columns]
        if cols: print(pre_df[cols].to_string(index=False))
    if not combo_df.empty:
        top10 = combo_df.sort_values('F1', ascending=False).head(10)
        top10_cols = [c for c in ['Subset', 'Strategy', 'N_Tools', 'Precision', 'Recall', 'F1', 'AUPRC']
                      if c in combo_df.columns]
        top10[top10_cols].to_csv(os.path.join(outpath, 'identification_best_combinations.tsv'), sep='\t', index=False)
        try:
            print(f"\n  Top 10 最佳组合 (F1):")
            print(top10[['Subset', 'Strategy', 'Precision', 'Recall', 'F1']].to_string(index=False))
            best = combo_df.loc[combo_df['F1'].astype(float).idxmax()]
            print(f"\n  最优 F1: {best['Subset']} ({best['Strategy']}) "
                  f"P={float(best['Precision']):.4f} R={float(best['Recall']):.4f} F1={float(best['F1']):.4f}")
        except Exception:
            pass
    return combo_df


# ── 核心处理流程 ────────────────────────────────────────────

def process_one_mode(labels_pd, y_true_all, engine, args, mode):
    """处理单个过滤模式，返回所有评估结果"""
    mode_label = FILTER_LABELS[mode]
    print(f"\n{'='*55}")
    print(f"  Filter Mode: {mode_label}")
    print(f"{'='*55}")
    if mode == 'raw':
        sub_dir = args.result_dir; suffix = args.suffix
    elif mode == 'filter':
        sub_dir = os.path.join(args.result_dir, 'uniprot_filter_output_filter')
        suffix = '.uniprot_filtered.id'
    else:
        sub_dir = os.path.join(args.result_dir, 'uniprot_filter_output_strict')
        suffix = '.uniprot_filtered.id'

    tool_keys = args.tools.split(',')
    if args.ensemble: tool_keys.append('all')

    all_overall, all_type_rows, cov_rows, strat_rows, entity_rows, entity_by_nsegs, fp_rows = [], [], [], [], [], [], []
    ds_store, y_preds_out = {}, {}

    for tool in tool_keys:
        fname = f'{args.prefix}.{tool}{suffix}' if tool != 'all' else f'{args.prefix}.all{suffix}'
        fpath = os.path.join(sub_dir, fname)
        if not os.path.exists(fpath):
            print(f"  [SKIP] 未找到: {fname}")
            continue

        base_label = TOOL_NAMES.get(tool, tool.upper())
        tool_label = base_label if mode == 'raw' else f"{base_label} ({mode_label})"

        metrics, y_pred = evaluate_tool(fpath, labels_pd, labels_pd, engine, tool_label)
        metrics['Mode'] = mode_label; all_overall.append(metrics)
        print(f"  [{tool_label:>28}] P={metrics['Precision']:.4f}  "
              f"R={metrics['Recall']:.4f}  F1={metrics['F1']:.4f}  N_Pred={metrics['N_Pred']}")

        ap, _, prec, rec = compute_auprc(y_true_all, y_pred.astype(float))
        ds_store[tool_label] = (prec, rec, ap)
        all_overall[-1]['AUPRC'] = round(ap, 4)
        y_preds_out[tool_label] = y_pred

        # 分层评估
        all_type_rows.extend(by_type_eval(labels_pd, y_pred, tool_label))
        cov_rows.extend(by_coverage_eval(labels_pd, y_pred, tool_label))
        fp_rows.append(extract_fp_breakdown(labels_pd, y_pred, tool_label))
        strat_rows.extend(by_scheme_stratified_eval(labels_pd, y_pred, tool_label))

        # 实体级评估
        ent_main, ent_nsegs = entity_level_eval(labels_pd, y_pred, tool_label)
        entity_rows.append(ent_main)
        entity_by_nsegs.extend(ent_nsegs)

    overall_df = pd.DataFrame(all_overall)
    if not overall_df.empty:
        overall_df = overall_df[overall_df['Precision'].notna()].drop_duplicates(subset='Tool')
    return {
        'overall': overall_df,
        'type_rows': all_type_rows,
        'cov_rows': cov_rows,
        'entity_rows': entity_rows,
        'entity_by_nsegs': entity_by_nsegs,
        'fp_rows': fp_rows,
        'strat_rows': strat_rows,
        'ds_store': ds_store,
        'y_preds': y_preds_out,
    }


# ── 主函数 ───────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='病毒鉴定评估脚本 (Polars加速+实体评估+FP细化 完整版)')
    parser.add_argument('--result-dir', default='step8_result/step3_identification_eval')
    parser.add_argument('--labels', required=True)
    parser.add_argument('--outdir', default='step8_result_analysis')
    parser.add_argument('--tools', default='blast,genomad,rdrpcatch,viralm,virbot,metabuli,viralverify,virhunter,virsorter2')
    parser.add_argument('--ensemble', action='store_true', default=True)
    parser.add_argument('--prefix', default='step3_identification_eval_virus')
    parser.add_argument('--suffix', default='.result.id')
    parser.add_argument('--filter-mode', default='raw', choices=['raw', 'filter', 'strict', 'all'])
    parser.add_argument('--virus-dir', default=None,
                        help='step1_eval_viruses/ 目录，用于排除短源病毒的阳性序列')
    parser.add_argument('--min-virus-length', type=int, default=1000,
                        help='排除源病毒基因组 < N bp 的阳性序列 (默认1000, 设为0禁用)')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()

    # ── 规范化输出子目录 (TSV 与对应图同目录) ──
    DIRS = {
        'overall':      os.path.join(args.outdir, 'overall'),       # PR曲线、ABCD四合一、汇总表、雷达图
        'modes':        os.path.join(args.outdir, 'modes'),         # 多模式对比: FP减少、散点、3x1柱状
        'overlap':      os.path.join(args.outdir, 'overlap'),       # 工具重叠分析
        'combinations': os.path.join(args.outdir, 'combinations'),  # 组合分析
        'coverage':     os.path.join(args.outdir, 'coverage'),      # 覆盖率召回
        'entity':       os.path.join(args.outdir, 'entity'),        # 实体评估 + FP拆解
    }
    for d in DIRS.values():
        os.makedirs(d, exist_ok=True)

    # ── 数据加载 (Polars 加速) ──
    engine = DataEngine()
    labels_pl = engine.load_labels(args.labels)
    # 排除短源病毒 (< min_virus_length bp) 的阳性序列
    if args.virus_dir and args.min_virus_length > 0:
        labels_pl = engine.filter_short_source_viruses(args.virus_dir, args.min_virus_length)
    labels_pd = labels_pl.to_pandas() if HAS_POLARS else labels_pl
    y_true_all = engine.get_y_true()
    modes_to_run = FILTER_MODES if args.filter_mode == 'all' else [args.filter_mode]

    # ── 逐模式处理 ──
    results = {}
    for mode in modes_to_run:
        results[mode] = process_one_mode(labels_pd, y_true_all, engine, args, mode)
        if results[mode]['overall'].empty:
            continue
        # 保存 TSV 到对应目录
        r = results[mode]
        r['overall'].to_csv(os.path.join(DIRS['overall'], f'identification_overall_{mode}.tsv'), sep='\t', index=False)
        pd.DataFrame(r['type_rows']).to_csv(os.path.join(DIRS['overall'], f'identification_by_type_{mode}.tsv'), sep='\t', index=False)
        pd.DataFrame(r['cov_rows']).to_csv(os.path.join(DIRS['coverage'], f'identification_by_coverage_{mode}.tsv'), sep='\t', index=False)
        pd.DataFrame(r['fp_rows']).to_csv(os.path.join(DIRS['entity'], f'fp_breakdown_{mode}.tsv'), sep='\t', index=False)
        pd.DataFrame(r['entity_rows']).to_csv(os.path.join(DIRS['entity'], f'entity_level_recall_{mode}.tsv'), sep='\t', index=False)
        if r['entity_by_nsegs']:
            pd.DataFrame(r['entity_by_nsegs']).to_csv(os.path.join(DIRS['entity'], f'entity_recall_by_nsegs_{mode}.tsv'), sep='\t', index=False)
        # MASTER scheme分层
        pd.DataFrame(r['strat_rows']).to_csv(os.path.join(DIRS['overall'], f'identification_stratified_{mode}.tsv'), sep='\t', index=False)

    if not results:
        print("未找到任何结果文件，退出。"); return

    # ── 跨模式聚合 ──
    all_overall = [r['overall'] for r in results.values() if not r['overall'].empty]
    if not all_overall:
        print("所有模式均无有效数据，退出。"); return

    raw_mode = 'raw' if 'raw' in results else modes_to_run[0]
    raw_overall = results[raw_mode]['overall']
    raw_ds = results[raw_mode]['ds_store']
    ds_stores = {m: results[m]['ds_store'] for m in modes_to_run}
    y_preds_per_mode = {m: results[m]['y_preds'] for m in modes_to_run}
    cov_df_all = pd.concat([pd.DataFrame(r['cov_rows']) for r in results.values() if r['cov_rows']], ignore_index=True)

    # ── 组合分析 (Raw 模式) ──
    combo_df_raw = pd.DataFrame()
    raw_yp = y_preds_per_mode.get(raw_mode, {})
    if raw_yp:
        tools_sorted = sorted([k for k in raw_yp.keys() if k != 'Ensemble'], key=lambda t: -np.sum(raw_yp[t]))
        combo_df_raw = run_combination_analysis(raw_yp, labels_pd, y_true_all, DIRS['combinations'], tools_sorted)

    # 三模式最优组合对比
    if len(modes_to_run) > 1:
        best_per_mode = []
        for mode in modes_to_run:
            yp = y_preds_per_mode.get(mode, {})
            if not yp: continue
            ts = sorted([k for k in yp.keys() if k != 'Ensemble'], key=lambda t: -np.sum(yp[t]))
            combo_df = run_combination_analysis(yp, labels_pd, y_true_all, DIRS['combinations'], ts)
            if not combo_df.empty:
                best = combo_df.loc[combo_df['F1'].astype(float).idxmax()]
                best_per_mode.append({'Mode': FILTER_LABELS[mode], **best})
        if best_per_mode:
            best_df = pd.DataFrame(best_per_mode)
            best_df.to_csv(os.path.join(DIRS['combinations'], 'identification_best_combo_by_mode.tsv'), sep='\t', index=False)

    # ── 制图阶段 ──
    if not args.no_plot:
        print(f"\n{'='*55}")
        print(f"  生成评估图表 ...")
        print(f"{'='*55}")

        # PR曲线 3×1: 三模式对比 (替代原单模式PR + 柱状图)
        if len(modes_to_run) > 1:
            plot_pr_curves_3x1(ds_stores, modes_to_run, os.path.join(DIRS['overall'], 'Fig_PR_Curves_3x1.png'))
        elif raw_ds:
            plot_pr_curves(raw_ds, os.path.join(DIRS['overall'], 'Fig_PR_Curves.png'))

        # Raw模式 ABCD 四合一概览 (仅作为初始分析概况)
        if raw_ds:
            cov_raw = cov_df_all[cov_df_all['Tool'].apply(lambda t: 'Filtered' not in t and 'Strict' not in t)] if len(cov_df_all) else cov_df_all
            plot_combined_figure(raw_ds, raw_overall, cov_raw,
                                 os.path.join(DIRS['overall'], 'Fig_Combined_ABCD_Raw.png'))

        # 汇总表图
        tbl_cols = [c for c in ['Tool', 'Precision', 'Recall', 'F1', 'AUPRC', 'MCC', 'Accuracy']
                    if c in raw_overall.columns]
        if tbl_cols:
            set_style(0.85)
            fig_tbl, ax_tbl = plt.subplots(figsize=(12, 4)); ax_tbl.axis('off')
            tbl_data = raw_overall[tbl_cols].round(4)
            tbl = ax_tbl.table(cellText=tbl_data.values, colLabels=tbl_data.columns,
                               cellLoc='center', loc='center')
            tbl.auto_set_font_size(False); tbl.set_fontsize(11); tbl.scale(1.2, 1.8)
            for j in range(len(tbl_cols)):
                tbl[0, j].set_text_props(fontweight='bold')
                tbl[0, j].set_facecolor('#4472C4')
                tbl[0, j].set_text_props(color='white')
            fig_tbl.savefig(os.path.join(DIRS['overall'], 'Fig_Metrics_Summary_Table.png'), bbox_inches='tight')
            plt.close(fig_tbl); print(f"  [overall] Fig_Metrics_Summary_Table.png")

        # 实体级评估 + FP拆解 → entity/ (3×1 合并图)
        if len(modes_to_run) > 1:
            all_fp_dfs = [pd.DataFrame(results[m]['fp_rows']) for m in modes_to_run]
            plot_fp_breakdown_dual_3x1(all_fp_dfs, modes_to_run, DIRS['entity'])
            all_ent_dfs = [pd.DataFrame(results[m]['entity_rows']) for m in modes_to_run]
            plot_entity_dumbbell_3x1(all_ent_dfs, modes_to_run, DIRS['entity'])
        else:
            r = results[modes_to_run[0]]
            fp_df = pd.DataFrame(r['fp_rows'])
            if not fp_df.empty:
                plot_fp_breakdown_dual(fp_df, DIRS['entity'], FILTER_LABELS[modes_to_run[0]])
            ent_df = pd.DataFrame(r['entity_rows'])
            if not ent_df.empty:
                plot_entity_dumbbell(ent_df, DIRS['entity'], FILTER_LABELS[modes_to_run[0]])

        # UpSet图 → overlap/ (每种模式单独 + 3×1 合并)
        for mode in modes_to_run:
            yp = y_preds_per_mode.get(mode, {})
            if yp:
                plot_overlap(yp, os.path.join(DIRS['overlap'], f'Fig_UpSet_{mode}.png'))
        if len(modes_to_run) > 1:
            plot_overlap_3x1(y_preds_per_mode, modes_to_run, os.path.join(DIRS['overlap'], 'Fig_UpSet_3x1.png'))

        # 覆盖率图 → coverage/
        if len(cov_df_all):
            plot_coverage_recall(cov_df_all, os.path.join(DIRS['coverage'], 'Fig_Recall_by_Coverage.png'))

        # MASTER scheme分层图 → overall/
        strat_all = pd.concat([pd.DataFrame(results[m]['strat_rows']) for m in modes_to_run if results[m].get('strat_rows')], ignore_index=True)
        if len(strat_all):
            for metric in ['Recall', 'Precision', 'F1', 'MCC', 'AUPRC']:
                if metric in strat_all.columns:
                    plot_scheme_stratified(strat_all, os.path.join(DIRS['overall'], f'Fig_Scheme_Stratified_{metric}.png'), metric)
            plot_scheme_radar_3x1(strat_all, os.path.join(DIRS['overall'], 'Fig_Scheme_Radar_3x1.png'))
            overall_by_scheme = strat_all[strat_all['Stratum'] == 'Overall']
            if len(overall_by_scheme):
                for metric in ['Recall', 'Precision', 'F1', 'MCC', 'AUPRC']:
                    if metric in overall_by_scheme.columns:
                        tbl = overall_by_scheme.pivot_table(index='Tool', columns='Scheme', values=metric, aggfunc='first')
                        tbl.to_csv(os.path.join(DIRS['overall'], f'table_scheme_{metric}_comparison.tsv'), sep='\t')

        # 组合分析图 → combinations/ (UpSet矩阵 + 边际收益 + 协同热力图)
        if not combo_df_raw.empty:
            plot_combinations(combo_df_raw, DIRS['combinations'])
            plot_marginal_benefit(combo_df_raw, DIRS['combinations'])
            plot_synergy_heatmap(combo_df_raw, DIRS['combinations'])

        # 多模式对比图 (仅在全模式时绘制)
        if len(modes_to_run) > 1:
            combined = pd.concat(all_overall, ignore_index=True)
            combined.to_csv(os.path.join(DIRS['modes'], 'identification_filter_comparison.tsv'), sep='\t', index=False)
            plot_fp_reduction(combined, os.path.join(DIRS['modes'], 'Fig_FP_Reduction.png'))
            plot_mode_scatter(combined, os.path.join(DIRS['modes'], 'Fig_Mode_Scatter.png'))
            plot_filter_comparison_3x1(combined, os.path.join(DIRS['modes'], 'Fig_Filter_Comparison_3x1.png'))
            plot_radar_charts(combined, DIRS['overall'])

    # ── 输出目录概览 ──
    print(f"\n{'='*55}")
    print(f"  输出目录结构 (TSV 与对应图同目录):")
    print(f"    {DIRS['overall']}/       — 总体评估 (PR、ABCD、柱状图、雷达图、汇总表)")
    print(f"    {DIRS['modes']}/         — 多模式对比 (FP减少、散点、3x1柱状)")
    print(f"    {DIRS['overlap']}/       — 工具重叠分析")
    print(f"    {DIRS['combinations']}/  — 组合分析 (热力图、Top20)")
    print(f"    {DIRS['coverage']}/      — 覆盖率召回曲线")
    print(f"    {DIRS['entity']}/        — 实体评估 + 假阳性来源拆解")
    print(f"{'='*55}")
    print(f"\n 评估完成!")


if __name__ == '__main__':
    main()
