#!/usr/bin/env python3
"""
宿主预测三工具评估：RVH vs PhaBOX2 vs C9_ICTV + Ensemble 集成

用法:
  python eval_host_prediction.py \
      --rvh step11-host/RVH_result/result.csv \
      --phabox step11-host/phabox2_output/final_prediction/cherry_prediction.tsv \
      --c9 step11-host/C9_ICTV_result/classification_result.tsv \
      --ensemble step11-host/ensemble_host_summary.tsv \
      --labels step3_identification_eval/sequence_labels_category.tsv \
      --outdir step11-host/evaluation/
"""

import argparse, os, re, sys
import pandas as pd
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

PALETTE = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F']

# ── 工具函数 ──────────────────────────────

def _extract_acc(s):
    """统一从 seq_id / contig_id 提取 accession token"""
    s = str(s).strip()
    m = re.search(r'(?:source=|src=)([^|]+)', s)
    if m: return m.group(1).strip()
    m = re.search(r'([A-Z]{1,4}_?\d{4,}\.\d{1,2})', s)
    if m: return m.group(1)
    return s

def _is_plant(val):
    v = str(val).strip().lower()
    return 'Plant' if v in ('plant', 'viridiplantae', 'streptophyta', 'algae') else 'Non-Plant'

# ── 数据加载 ──────────────────────────────

def load_labels(path):
    df = pd.read_csv(path, sep='\t')
    df['Coverage'] = df['seq_id'].str.extract(r'cov(\d+)')[0].astype(float)
    print(f"[Labels] {len(df)} seqs | Pos={sum(df['label']=='positive')}")
    return df

def load_rvh(path):
    df = pd.read_csv(path)
    df['contig_id'] = df.iloc[:, 0].apply(str).str.strip()
    df['Host_Pred'] = df['pred|L1'].fillna('Unknown')
    df['Is_Plant'] = df['Host_Pred'].apply(_is_plant)
    df['_key'] = df['contig_id'].apply(_extract_acc)  # 用于 C9/Ens 交叉比对
    return df.set_index('contig_id')  # RVH 用完整 seq_id 做 index

def load_phabox(path):
    df = pd.read_csv(path, sep='\t').rename(columns={'Accession': 'contig_id'})
    df['Host_Pred'] = df['Host'].fillna('Unknown')
    df['Is_Plant'] = df['Host_Pred'].apply(_is_plant)
    df['_key'] = df['contig_id'].apply(_extract_acc)
    return df.set_index('contig_id')

def load_c9(path):
    df = pd.read_csv(path, sep='\t')
    df['Host_Pred'] = df['Predicted_Host'].fillna('Unknown')
    df['Is_Plant'] = df['Host_Pred'].apply(_is_plant)
    df['_key'] = df['contig_id'].apply(_extract_acc)
    # C9 多行共享一个 accession → 取多数投票
    grp = df.groupby('_key').agg(
        Plant_votes=('Is_Plant', lambda x: sum(x == 'Plant')),
        Total_votes=('Is_Plant', 'count'),
        Host_Pred=('Host_Pred', 'first')
    ).reset_index()
    grp['Is_Plant'] = grp.apply(lambda r: 'Plant' if r['Plant_votes'] > r['Total_votes'] / 2 else 'Non-Plant', axis=1)
    print(f"  [load_c9] {len(df)} rows → {len(grp)} unique keys (majority vote), Plant={sum(grp['Is_Plant']=='Plant')}")
    return grp.set_index('_key')

def load_ensemble(path):
    df = pd.read_csv(path, sep='\t')
    df['Host_Pred'] = df['Final_Host'].fillna('Unknown')
    df['Is_Plant'] = df['Host_Pred'].apply(_is_plant)
    df['_key'] = df['contig_id'].apply(_extract_acc)
    grp = df.groupby('_key').agg(
        Plant_votes=('Is_Plant', lambda x: sum(x == 'Plant')),
        Total_votes=('Is_Plant', 'count'),
    ).reset_index()
    grp['Is_Plant'] = grp.apply(lambda r: 'Plant' if r['Plant_votes'] > r['Total_votes'] / 2 else 'Non-Plant', axis=1)
    print(f"  [load_ens] {len(df)} rows → {len(grp)} unique keys (majority vote), Plant={sum(grp['Is_Plant']=='Plant')}")
    return grp.set_index('_key')

# ── 评估 ──────────────────────────────────

def evaluate(df_pred, labels_df, tool_name):
    labels_df['_key'] = labels_df['seq_id'].apply(_extract_acc)
    if '_key' in df_pred.columns:
        pred_keys = set(df_pred['_key'])
        df_pred_flat = df_pred[['_key', 'Is_Plant']].drop_duplicates(subset='_key')
    else:
        pred_keys = set(df_pred.index)
        df_pred_flat = df_pred.reset_index()[['_key', 'Is_Plant']]
    common_keys = set(labels_df['_key']) & pred_keys
    if len(common_keys) == 0:
        print(f'  [{tool_name}] No matching keys (pred={len(pred_keys)}, label={labels_df["_key"].nunique()})')
        return {}
    lab_grp = labels_df[labels_df['_key'].isin(common_keys)].groupby('_key').agg(n_pos=('label', lambda x: sum(x == "positive")), n_total=('label', 'count')).reset_index()
    merged = lab_grp.merge(df_pred_flat[df_pred_flat['_key'].isin(common_keys)], on='_key', how='inner')
    if len(merged) == 0:
        print(f'  [{tool_name}] No entries after merge')
        return {}
    y_true = (merged['n_pos'] > 0).astype(int).values
    y_pred = (merged['Is_Plant'] == 'Plant').astype(int).values
    TP = int(sum((y_true == 1) & (y_pred == 1)))
    FP = int(sum((y_true == 0) & (y_pred == 1)))
    FN = int(sum((y_true == 1) & (y_pred == 0)))
    TN = int(sum((y_true == 0) & (y_pred == 0)))
    precision = TP / max(1, TP + FP)
    recall = TP / max(1, TP + FN)
    f1 = 2 * precision * recall / max(1e-10, precision + recall)
    accuracy = (TP + TN) / max(1, len(y_true))
    print(f'  [{tool_name:>12}] Keys={len(merged):>4} | P={precision:.4f} R={recall:.4f} F1={f1:.4f} Acc={accuracy:.4f}')
    return {'Tool': tool_name, 'Keys': len(merged), 'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN, 'Precision': round(precision, 4), 'Recall': round(recall, 4), 'F1': round(f1, 4), 'Accuracy': round(accuracy, 4), 'Plant_Predictions': int(sum(y_pred)), 'Total_Plant': int(sum(y_true))}

def coverage_recall(df_pred, labels_df, tool_name):
    pos = labels_df[labels_df['label'] == 'positive'].copy()
    pos['_key'] = pos['seq_id'].apply(_extract_acc)
    if '_key' in df_pred.columns:
        pred_keys = set(df_pred['_key'])
        df_pred_flat = df_pred[['_key', 'Is_Plant']].drop_duplicates(subset='_key')
    else:
        pred_keys = set(df_pred.index)
        df_pred_flat = df_pred.reset_index()[['_key', 'Is_Plant']]
    common = set(pos['_key']) & pred_keys
    if len(common) == 0:
        return []
    rows = []
    for cov in sorted(pos['Coverage'].dropna().unique()):
        cov_pos = pos[(pos['Coverage'] == cov) & (pos['_key'].isin(common))]
        cov_keys = set(cov_pos['_key'])
        pred_sub = df_pred_flat[df_pred_flat['_key'].isin(cov_keys)]
        n = len(cov_keys)
        if n == 0: continue
        plant_hits = sum(pred_sub['Is_Plant'] == 'Plant') if len(pred_sub) > 0 else 0
        rows.append({'Coverage': f'{int(cov)}%', 'N_Seqs': n, 'Detected': int(plant_hits), 'Recall': round(int(plant_hits) / max(1, n), 4), 'Tool': tool_name})
    return rows

def plot_comparison(metrics_df, outpath):
    sns.set_theme(style='whitegrid', font_scale=1.1)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric in zip(axes, ['Precision', 'Recall', 'F1']):
        bars = ax.bar(metrics_df['Tool'], metrics_df[metric].astype(float),
                      color=['#4C72B0', '#55A868', '#C44E52', '#7E6148'], edgecolor='white')
        ax.set_title(metric, fontsize=14); ax.set_ylim(0, 1.15)
        for bar, val in zip(bars, metrics_df[metric]):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f'{val:.3f}', ha='center')
        ax.tick_params(axis='x', rotation=30)
    fig.suptitle('Host Prediction: Plant vs Non-Plant', fontweight='bold')
    plt.tight_layout(); fig.savefig(outpath, dpi=300); plt.close(fig)

def plot_confusion_matrices(metrics_list, outpath):
    fig, axes = plt.subplots(1, len(metrics_list), figsize=(4 * len(metrics_list), 4))
    if len(metrics_list) == 1: axes = [axes]
    for ax, m in zip(axes, metrics_list):
        cm = np.array([[m['TP'], m['FN']], [m['FP'], m['TN']]])
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax,
                    xticklabels=['Plant', 'Non-Plant'], yticklabels=['Plant', 'Non-Plant'])
        ax.set_title(m['Tool'], fontweight='bold')
        ax.set_xlabel('Predicted'); ax.set_ylabel('True')
    plt.suptitle('Confusion Matrices', fontweight='bold')
    plt.tight_layout(); fig.savefig(outpath, dpi=300); plt.close(fig)

def plot_coverage(cov_df, outpath):
    if cov_df.empty: return
    sns.set_theme(style='whitegrid'); fig, ax = plt.subplots(figsize=(9, 5))
    for i, tool in enumerate(cov_df['Tool'].unique()):
        sub = cov_df[cov_df['Tool'] == tool]
        ax.plot(sub['Coverage'], sub['Recall'].astype(float), marker='o', linewidth=2,
                color=PALETTE[i % len(PALETTE)], label=tool)
    ax.set_xlabel('Coverage'); ax.set_ylabel('Recall')
    ax.set_title('Plant Host Recall by Coverage', fontweight='bold')
    ax.legend(); ax.set_ylim(0, 1.1); ax.grid(alpha=0.3)
    plt.tight_layout(); fig.savefig(outpath, dpi=300); plt.close(fig)

def plot_radar(metrics_df, outpath):
    """雷达图: 4工具 × Precision/Recall/F1/Accuracy"""
    from math import pi
    categories = ['Precision', 'Recall', 'F1', 'Accuracy']
    N = len(categories)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]  # 闭合

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    for i, tool in enumerate(metrics_df['Tool']):
        values = [float(metrics_df[metrics_df['Tool'] == tool][c].iloc[0]) for c in categories]
        values += values[:1]
        ax.fill(angles, values, alpha=0.15, color=PALETTE[i % len(PALETTE)])
        ax.plot(angles, values, 'o-', linewidth=2, color=PALETTE[i % len(PALETTE)], label=tool)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=12)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=8, color='grey')
    ax.set_title('Host Prediction: Multi-Metric Radar', fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    plt.tight_layout(); fig.savefig(outpath, dpi=300); plt.close(fig)

# ── 主函数 ────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='宿主预测三工具评估')
    parser.add_argument('--rvh', help='RVH result.csv')
    parser.add_argument('--phabox', help='PhaBOX2 cherry_prediction.tsv')
    parser.add_argument('--c9', help='C9_ICTV classification_result.tsv')
    parser.add_argument('--ensemble', help='ensemble_host_summary.tsv')
    parser.add_argument('--labels', required=True)
    parser.add_argument('--outdir', default='host_evaluation')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    labels = load_labels(args.labels)
    all_metrics, all_cov = [], []

    for path, loader, name in [
        (args.rvh,      load_rvh,      'RVH'),
        (args.phabox,   load_phabox,   'PhaBOX2'),
        (args.c9,       load_c9,       'C9_ICTV'),
        (args.ensemble, load_ensemble, 'Ensemble'),
    ]:
        if not path or not os.path.exists(path):
            print(f"  [SKIP] {path} not found")
            continue
        df_pred = loader(path)
        m = evaluate(df_pred, labels, name)
        if m:
            all_metrics.append(m)
            all_cov.extend(coverage_recall(df_pred, labels, name))

    metrics_df = pd.DataFrame(all_metrics)
    if metrics_df.empty:
        print("No metrics generated. Check input paths.")
        sys.exit(1)

    m_path = os.path.join(args.outdir, 'host_prediction_metrics.tsv')
    metrics_df.to_csv(m_path, sep='\t', index=False)
    print(f"\n✅ {m_path}")
    print(metrics_df.to_string(index=False))

    cov_df = pd.DataFrame(all_cov)
    if len(cov_df):
        cov_df.to_csv(os.path.join(args.outdir, 'host_coverage_recall.tsv'), sep='\t', index=False)

    plot_comparison(metrics_df, os.path.join(args.outdir, 'Fig_Host_Comparison.png'))
    plot_confusion_matrices(all_metrics, os.path.join(args.outdir, 'Fig_Host_Confusion.png'))
    plot_radar(metrics_df, os.path.join(args.outdir, 'Fig_Host_Radar.png'))
    if len(cov_df):
        plot_coverage(cov_df, os.path.join(args.outdir, 'Fig_Host_Coverage_Recall.png'))

    print(f"\n🎉 Done: {args.outdir}")


if __name__ == '__main__':
    main()
