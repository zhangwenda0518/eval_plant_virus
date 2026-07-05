#!/usr/bin/env python3
"""
去冗余聚类评估 + 可视化
金标准：片段 ID 中的物种名（如 NC_002030_mut0pct_len50pct_f1 → NC_002030）
评估指标：ARI, AMI, Homogeneity, Completeness, V-measure, Purity, NMI
分层分析：按突变率、按长度比例

用法:
  python eval_dedup_clustering.py \
      --input step2_dedup_fragments.fasta \
      --cluster-dir step3_dedup_cluster/ \
      --outdir step4_dedup_eval/ --min-id 0.90
"""

import argparse, os, re, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from sklearn.metrics import (adjusted_rand_score, adjusted_mutual_info_score,
                              homogeneity_score, completeness_score,
                              v_measure_score, normalized_mutual_info_score)

# ── 配置 ──
MY_PALETTE = {'MMseqs2': '#4C72B0', 'CD-HIT': '#55A868', 'VCLUST': '#DD8452',
              'Linclust': '#A0522D'}

MUT_ORDER = ['0%', '5%', 'host']
LEN_ORDER = ['20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%']


def parse_gold_from_fasta(fasta_path):
    """从 FASTA ID 解析金标准物种标签和突变/长度属性
    病毒: NC_002030_mut0pct_len50pct_f42_pos... → species=NC_002030
    宿主: HOST_host_f1_len50pct_pos...          → species=HOST_host_f1 (unique per region)
    """
    labels = {}
    attrs = {}
    for line in open(fasta_path):
        if not line.startswith('>'):
            continue
        seq_id = line[1:].strip().split()[0]
        parts = seq_id.split('_')

        if parts[0] == 'HOST':
            # 宿主片段：每个 host_f{N} 是一个独立的 gold singleton
            # ID: HOST_host_f1_len50pct_pos0-3000_sw0-2500
            species = '_'.join(parts[:2])  # HOST_host_f1
            labels[seq_id] = species
            len_str = "?"
            for p in parts:
                if p.startswith('len'):
                    len_str = p.replace('len', '').replace('pct', '%')
            attrs[seq_id] = {'mutation': 'host', 'length_fraction': len_str}
        else:
            species = parts[0]
            labels[seq_id] = species
            mut_str = "0%"; len_str = "?"
            for p in parts:
                if p.startswith('mut'):
                    mut_str = p.replace('mut', '').replace('pct', '%')
                if p.startswith('len'):
                    len_str = p.replace('len', '').replace('pct', '%')
            attrs[seq_id] = {'mutation': mut_str, 'length_fraction': len_str}

    n_species = len(set(labels.values()))
    n_host = sum(1 for v in labels.values() if v.startswith('HOST'))
    print(f"[gold] {len(labels)} seqs, {n_species} species ({n_host} host singletons)")
    return labels, attrs


def parse_mmseqs_clusters(cluster_tsv):
    """MMseqs2 easy-cluster: rep \t member"""
    clusters = defaultdict(set)
    for line in open(cluster_tsv):
        rep, member = line.strip().split('\t')
        clusters[rep].add(rep)
        clusters[rep].add(member)
    return dict(clusters)


def parse_cdhit_clusters(clstr_file):
    """CD-HIT .clstr 格式"""
    clusters = defaultdict(set)
    current_cluster = None
    for line in open(clstr_file):
        if line.startswith('>'):
            current_cluster = line.strip().split()[1]
        elif current_cluster:
            seq_id = line.split('>')[1].split('...')[0].strip()
            clusters[current_cluster].add(seq_id)
    return dict(clusters)


def parse_vclust_clusters(uc_file):
    """VSEARCH .uc 格式"""
    clusters = defaultdict(set)
    for line in open(uc_file):
        parts = line.strip().split('\t')
        if parts[0] == 'C':
            cluster_id = parts[1]
            seq_id = parts[8]
            clusters[cluster_id].add(seq_id)
        elif parts[0] == 'H':
            cluster_id = parts[1]
            seq_id = parts[8]
            clusters[cluster_id].add(seq_id)
    return dict(clusters)


def compute_metrics(true_labels, pred_clusters):
    """ARI, AMI, Homogeneity, Completeness, V-measure, Purity"""
    all_ids = list(true_labels.keys())
    y_true = [true_labels[sid] for sid in all_ids]
    pred_map = {}
    for cluster_name, members in pred_clusters.items():
        for m in members:
            if m in pred_map:
                pred_map[m] = cluster_name  # 取首次出现
            else:
                pred_map[m] = cluster_name
    y_pred = [pred_map.get(sid, f'__unassigned_{hash(sid)}') for sid in all_ids]

    return {
        'ARI': adjusted_rand_score(y_true, y_pred),
        'AMI': adjusted_mutual_info_score(y_true, y_pred),
        'NMI': normalized_mutual_info_score(y_true, y_pred),
        'Homogeneity': homogeneity_score(y_true, y_pred),
        'Completeness': completeness_score(y_true, y_pred),
        'V_measure': v_measure_score(y_true, y_pred),
        'n_seqs': len(all_ids),
        'n_species': len(set(y_true)),
        'n_clusters': len(pred_clusters),
    }


def stratified_metrics(true_labels, attrs, pred_clusters, stratify_by):
    """按突变率或长度比例分层计算指标"""
    results = []
    strata = set(a[stratify_by] for a in attrs.values())
    for val in sorted(strata):
        subset_ids = {sid for sid, a in attrs.items() if a[stratify_by] == val}
        sub_true = {sid: true_labels[sid] for sid in subset_ids if sid in true_labels}
        sub_pred = {}
        for cname, members in pred_clusters.items():
            filtered = members & set(sub_true.keys())
            if filtered:
                sub_pred[cname] = filtered
        if not sub_true or not sub_pred:
            continue
        m = compute_metrics(sub_true, sub_pred)
        m[stratify_by] = val
        results.append(m)
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="去冗余聚类评估")
    parser.add_argument("--input", required=True, help="模拟片段 FASTA")
    parser.add_argument("--cluster-dir", required=True, help="聚类结果目录")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-id", type=float, default=0.90)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. 金标准
    true_labels, attrs = parse_gold_from_fasta(args.input)

    # 2. 解析各工具聚类结果
    tools = {}
    # MMseqs2
    mmseqs_files = glob.glob(os.path.join(args.cluster_dir, "mmseqs2", "*_cluster.tsv"))
    if mmseqs_files:
        tools['MMseqs2'] = parse_mmseqs_clusters(mmseqs_files[0])
    # CD-HIT
    cdhit_files = glob.glob(os.path.join(args.cluster_dir, "cdhit", "*.clstr"))
    if cdhit_files:
        tools['CD-HIT'] = parse_cdhit_clusters(cdhit_files[0])
    # VCLUST
    vclust_files = glob.glob(os.path.join(args.cluster_dir, "vclust", "*.uc"))
    if vclust_files:
        tools['VCLUST'] = parse_vclust_clusters(vclust_files[0])
    # Linclust
    lc_files = glob.glob(os.path.join(args.cluster_dir, "linclust", "*_cluster.tsv"))
    if lc_files:
        tools['Linclust'] = parse_mmseqs_clusters(lc_files[0])

    if not tools:
        print("[ERROR] 未找到任何聚类结果！")
        return

    print(f"[tools] Loaded: {list(tools.keys())}")

    # 3. 整体评估
    overall_rows = []
    stratified_mut = []
    stratified_len = []
    for name, clusters in tools.items():
        m = compute_metrics(true_labels, clusters)
        m['Tool'] = name
        overall_rows.append(m)
        # 按突变率分层
        mut_df = stratified_metrics(true_labels, attrs, clusters, 'mutation')
        mut_df['Tool'] = name
        stratified_mut.append(mut_df)
        # 按长度分层
        len_df = stratified_metrics(true_labels, attrs, clusters, 'length_fraction')
        len_df['Tool'] = name
        stratified_len.append(len_df)

    overall = pd.DataFrame(overall_rows)
    print("\n  🏆 Overall Performance:")
    for _, r in overall.iterrows():
        print(f"    {r['Tool']:12s}  ARI={r['ARI']:.4f}  AMI={r['AMI']:.4f}  "
              f"V={r['V_measure']:.4f}  NMI={r['NMI']:.4f}  n_clust={r['n_clusters']}")

    overall.to_csv(os.path.join(args.outdir, "dedup_overall.tsv"), sep='\t', index=False)

    # 4. 分层表
    mut_all = pd.concat(stratified_mut) if stratified_mut else pd.DataFrame()
    len_all = pd.concat(stratified_len) if stratified_len else pd.DataFrame()
    if not mut_all.empty:
        mut_all.to_csv(os.path.join(args.outdir, "dedup_by_mutation.tsv"), sep='\t', index=False)
    if not len_all.empty:
        len_all.to_csv(os.path.join(args.outdir, "dedup_by_length.tsv"), sep='\t', index=False)

    # 5. 可视化
    sns.set_theme(style='ticks', font_scale=1.05)
    for metric, ylabel in [('ARI', 'Adjusted Rand Index'), ('V_measure', 'V-measure')]:
        # 整体条形图
        fig, ax = plt.subplots(figsize=(6, 5))
        order = overall.sort_values(metric, ascending=False)['Tool']
        sns.barplot(data=overall, x='Tool', y=metric, order=order, palette=MY_PALETTE, ax=ax)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} by Tool (id={args.min_id})", fontweight='bold')
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"bar_{metric}.png"), dpi=200, bbox_inches='tight')
        plt.close()

        # 按突变率分层
        if not mut_all.empty:
            fig, ax = plt.subplots(figsize=(8, 5))
            sns.boxplot(data=mut_all, x='mutation', y=metric, hue='Tool',
                        order=MUT_ORDER, hue_order=list(tools.keys()),
                        palette=MY_PALETTE, ax=ax)
            ax.set_ylim(0, 1.02)
            ax.set_ylabel(ylabel)
            ax.set_title(f"{ylabel} by Mutation Rate", fontweight='bold')
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, f"box_{metric}_by_mutation.png"), dpi=200, bbox_inches='tight')
            plt.close()

        # 按长度比例分层
        if not len_all.empty:
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=len_all, x='length_fraction', y=metric, hue='Tool',
                        hue_order=list(tools.keys()), palette=MY_PALETTE,
                        markers=True, ax=ax)
            ax.set_ylim(0, 1.02)
            ax.set_ylabel(ylabel)
            ax.set_xlabel("Length Fraction")
            ax.set_title(f"{ylabel} by Fragment Length", fontweight='bold')
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, f"line_{metric}_by_length.png"), dpi=200, bbox_inches='tight')
            plt.close()

    print(f"\n✅ Done. Output: {args.outdir}")
    print(f"  dedup_overall.tsv  dedup_by_mutation.tsv  dedup_by_length.tsv")


if __name__ == "__main__":
    main()
