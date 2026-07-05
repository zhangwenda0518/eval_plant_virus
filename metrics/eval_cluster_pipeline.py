#!/usr/bin/env python3
"""
run_cluster_pipeline.py 输出评估：ARI/AMI + 分层 + 可视化

用法:
  python eval_cluster_pipeline.py \
      --pipeline-dir Clustering_Pipeline_Out \
      --outdir cluster_eval/ --min-id 0.95
"""

import argparse, os, re, glob
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import (adjusted_rand_score, adjusted_mutual_info_score,
    homogeneity_score, completeness_score, v_measure_score, normalized_mutual_info_score)

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

MY_PALETTE = {'mmseqs':'#55A868','vclust':'#DD8452','blast':'#4C72B0','drep':'#E69F00',
              'cdhit':'#8B5CF6','sumaclust':'#F39B7F','hdbscan':'#C44E52'}


def parse_gold(fasta_dir_or_file):
    """从 FASTA ID 提取 gold species + 属性"""
    labels, attrs = {}, {}
    files = list(glob.glob(os.path.join(fasta_dir_or_file, "*.fasta"))) if os.path.isdir(fasta_dir_or_file) \
            else [fasta_dir_or_file]
    for f in files:
        for rec in SeqIO.parse(f, "fasta"):
            sid = rec.id.split()[0]
            parts = sid.split('_')
            # NC_014509.2 类型的 accession 按下划线 split 后 parts[0] 只是前缀
            if len(parts) >= 2 and parts[0] in ('NC','AC','NG','NT','NW','NZ') and '.' in parts[1]:
                labels[sid] = parts[0] + '_' + parts[1]
            else:
                labels[sid] = parts[0]
            mut, lf = '?', '?'
            for p in parts:
                if p.startswith('mut'): mut = p.replace('mut', '').replace('pct', '%')
                if p.startswith('len'): lf = p.replace('len', '').replace('pct', '%')
            attrs[sid] = {'mutation': mut, 'length_fraction': lf}
    return labels, attrs


def parse_aniclust_clusters(tsv):
    """BLAST/dRep clusters.tsv: rep\\tmember1,member2,..."""
    clusters = defaultdict(set)
    for line in open(tsv):
        if line.startswith('representative'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            rep = parts[0]
            clusters[rep].add(rep)
            for m in parts[1].split(','):
                m = m.strip()
                if m:
                    clusters[rep].add(m)
    return dict(clusters)


def parse_pipeline_clusters(pipeline_dir, tool):
    """从 Split_Fastas 读取聚类结果"""
    sf_dir = os.path.join(pipeline_dir, f"{tool}_results", "Split_Fastas")
    clusters = defaultdict(set)
    for f in sorted(glob.glob(os.path.join(sf_dir, "*.all.fasta"))):
        cid = os.path.basename(f).replace('.all.fasta', '')
        for rec in SeqIO.parse(f, "fasta"):
            clusters[cid].add(rec.id.split()[0])
    return dict(clusters)


def compute_metrics(true_labels, pred_clusters):
    all_ids = list(true_labels.keys())
    y_true = [true_labels[i] for i in all_ids]
    pmap = {}
    for cn, ms in pred_clusters.items():
        for m in ms: pmap[m] = cn
    y_pred = [pmap.get(i, f'__u_{hash(i)}') for i in all_ids]
    # Purity: sum(max_count_per_cluster) / n (去重防序列跨簇重复)
    seen = set()
    purity_sum = 0
    for cn, ms in pred_clusters.items():
        c_true = [true_labels[m] for m in ms if m in true_labels and m not in seen]
        for m in ms:
            seen.add(m)
        if c_true:
            purity_sum += max(pd.Series(c_true).value_counts())
    purity = purity_sum / len(all_ids) if all_ids else 0

    return {
        'ARI': adjusted_rand_score(y_true, y_pred),
        'AMI': adjusted_mutual_info_score(y_true, y_pred),
        'NMI': normalized_mutual_info_score(y_true, y_pred),
        'Homogeneity': homogeneity_score(y_true, y_pred),
        'Completeness': completeness_score(y_true, y_pred),
        'V_measure': v_measure_score(y_true, y_pred),
        'Purity': round(purity, 4),
        'n_seqs': len(all_ids), 'n_species': len(set(y_true)),
        'n_clusters': len(pred_clusters),
    }


def stratified_metrics(true_labels, attrs, clusters, stratify_by):
    results = []
    strata = set(a[stratify_by] for a in attrs.values())
    for val in sorted(strata):
        sids = {i for i, a in attrs.items() if a[stratify_by] == val}
        sub_true = {i: true_labels[i] for i in sids if i in true_labels}
        sub_pred = {}; sset = set(sub_true.keys())
        for cn, ms in clusters.items():
            fm = ms & sset
            if fm: sub_pred[cn] = fm
        if sub_true and sub_pred:
            m = compute_metrics(sub_true, sub_pred)
            m[stratify_by] = val; results.append(m)
    return pd.DataFrame(results)


def plot(overall, out_dir, min_id):
    sns.set_theme(style='ticks', font_scale=1.1)
    for metric, yl in [('ARI', 'Adjusted Rand Index'), ('V_measure', 'V-measure')]:
        fig, ax = plt.subplots(figsize=(5, 4.5))
        order = overall.sort_values(metric, ascending=False)['Tool']
        sns.barplot(data=overall, x='Tool', y=metric, order=order, hue='Tool',
                    palette=MY_PALETTE, legend=False, ax=ax)
        ax.set_ylim(0, 1.02); ax.set_ylabel(yl)
        ax.set_title(f"{yl} (id={min_id})", fontweight='bold')
        for p in ax.patches:
            ax.annotate(f'{p.get_height():.3f}', (p.get_x()+p.get_width()/2, p.get_height()),
                       ha='center', va='bottom', fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"{metric}.png"), dpi=200, bbox_inches='tight')
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="聚类管道评估")
    parser.add_argument("--pipeline-dir", required=True, help="run_cluster_pipeline.py 输出目录")
    parser.add_argument("--input-fasta", default=None, help="原始模拟片段 FASTA (可选，默认从 pipeline out 读取)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-id", type=float, default=0.95, help="(仅标注用)")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Gold labels
    fasta_src = args.input_fasta or os.path.join(args.pipeline_dir, "mmseqs_results", "all.cluster.ref.fasta")
    true_labels, attrs = parse_gold(fasta_src)
    print(f"[gold] {len(true_labels)} seqs, {len(set(true_labels.values()))} species")

    # Parse each tool's clusters
    tools = {}
    for tool in ['mmseqs', 'vclust', 'blast', 'drep', 'cdhit', 'sumaclust', 'hdbscan']:
        td = os.path.join(args.pipeline_dir, f"{tool}_results")
        if os.path.isdir(td) and os.path.isdir(os.path.join(td, "Split_Fastas")):
            tools[tool] = parse_pipeline_clusters(args.pipeline_dir, tool)
        # BLAST/dRep/aniclust 输出 clusters.tsv (rep\tm1,m2,...)
        elif os.path.isdir(td):
            tsv = os.path.join(td, "clusters.tsv")
            if os.path.exists(tsv):
                tools[tool] = parse_aniclust_clusters(tsv)
        if tool in tools:
            print(f"  [{tool}] {len(tools[tool])} clusters")

    if not tools:
        print("[ERROR] No clustering results"); return

    # Overall
    overall_rows = []
    mut_rows, len_rows = [], []
    for name, clusters in tools.items():
        m = compute_metrics(true_labels, clusters)
        m['Tool'] = name; overall_rows.append(m)
        mdf = stratified_metrics(true_labels, attrs, clusters, 'mutation')
        if not mdf.empty: mdf['Tool'] = name; mut_rows.append(mdf)
        ldf = stratified_metrics(true_labels, attrs, clusters, 'length_fraction')
        if not ldf.empty: ldf['Tool'] = name; len_rows.append(ldf)

    overall = pd.DataFrame(overall_rows)
    print(f"\n  🏆 Overall (id≈{args.min_id}):")
    for _, r in overall.iterrows():
        print(f"    {r['Tool']:10s}  ARI={r['ARI']:.4f}  AMI={r['AMI']:.4f}  V={r['V_measure']:.4f}  "
              f"H={r['Homogeneity']:.4f}  C={r['Completeness']:.4f}  Purity={r['Purity']:.4f}  n_clust={r['n_clusters']}")

    overall.to_csv(os.path.join(args.outdir, "eval_overall.tsv"), sep='\t', index=False)
    mut_all = pd.concat(mut_rows) if mut_rows else pd.DataFrame()
    len_all = pd.concat(len_rows) if len_rows else pd.DataFrame()
    if not mut_all.empty: mut_all.to_csv(os.path.join(args.outdir, "eval_by_mutation.tsv"), sep='\t', index=False)
    if not len_all.empty: len_all.to_csv(os.path.join(args.outdir, "eval_by_length.tsv"), sep='\t', index=False)

    plot(overall, args.outdir, args.min_id)

    # 资源评估
    res_rows = []
    for tool in ['mmseqs', 'vclust']:
        for f in glob.glob(os.path.join(args.pipeline_dir, f"**/{tool}*.time.mem.log"), recursive=True):
            try:
                content = open(f).read()
                t = re.search(r'Time:([\d.]+)', content)
                m = re.search(r'Memory:(\d+)', content)
                if t:
                    res_rows.append({'Tool': tool, 'Time_s': float(t.group(1)),
                                     'Mem_MB': round(int(m.group(1))/1024, 1) if m else 0,
                                     'Step': os.path.basename(f).replace('.time.mem.log','')})
            except Exception:
                pass
    if res_rows:
        res_df = pd.DataFrame(res_rows)
        res_df.to_csv(os.path.join(args.outdir, "resource.tsv"), sep='\t', index=False)
        # 打印
        grp = res_df.groupby('Tool')[['Time_s', 'Mem_MB']].sum()
        print(f"\n  📊 Resources:")
        for tool in grp.index:
            r = grp.loc[tool]
            print(f"    {tool:10s}  Total Time={r['Time_s']:.1f}s  Peak Mem={r['Mem_MB']:.0f}MB")
        # 图
        fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
        for ax, col, yl in zip(axes, ['Time_s', 'Mem_MB'], ['Total Time (s)', 'Peak Memory (MB)']):
            sns.barplot(data=grp.reset_index(), x='Tool', y=col, hue='Tool',
                       palette=MY_PALETTE, legend=False, ax=ax)
            ax.set_ylabel(yl)
            for p in ax.patches:
                ax.annotate(f'{p.get_height():.0f}', (p.get_x()+p.get_width()/2, p.get_height()),
                           ha='center', va='bottom', fontsize=9)
        fig.suptitle('Resource Comparison', fontweight='bold')
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, "resource.png"), dpi=200, bbox_inches='tight')
        plt.close()
        print(f"   📊 resource.tsv, resource.png")

    print(f"\n✅ Done. Output: {args.outdir}")


if __name__ == "__main__":
    main()
