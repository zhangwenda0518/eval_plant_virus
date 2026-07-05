#!/usr/bin/env python3
"""
构建按病毒科排除的未知病毒鉴定评估数据集 (Leave-One-Family-Out, 参考 VirHunter 2022)

策略:
  1. 从 selected_viruses.tsv 获取每个病毒的科 (Family) 信息
  2. 对每个包含 ≥2 个病毒的科, 将其全部病毒标记为该轮的 "未知病毒"
  3. 按覆盖率梯度截取片段作为正样本
  4. 输出每科的独立评估文件, 以及合并的 "全未知" 文件

用法:
  python prep_build_family_leaveout_eval.py \
      --virus-dir step1_eval_viruses/ \
      --virus-meta step1_eval_viruses/selected_viruses.tsv \
      --coverage-levels 100 80 60 40 \
      --n-per-coverage 2 \
      --host-genome ningxia.genome.fasta \
      --conserved-prots step3_conserved_traps/ \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --n-neg-A 200 --n-neg-B 150 --n-neg-C 100 \
      --outdir step3_family_leaveout_eval/ --seed 42
"""

import argparse, os, sys, random
from pathlib import Path
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np


def load_virus_metadata(meta_path):
    """加载 selected_viruses.tsv, 返回 {accession: {family, genus, species, length}}"""
    meta = pd.read_csv(meta_path, sep='\t')
    info = {}
    for _, row in meta.iterrows():
        acc = row['accession']
        info[acc] = {
            'family': str(row['family']) if pd.notna(row.get('family')) else 'Unclassified',
            'genus': str(row['genus']) if pd.notna(row.get('genus')) else 'Unknown',
            'species': str(row['species']) if pd.notna(row.get('species')) else 'Unknown',
            'length': int(row['length']) if pd.notna(row.get('length')) else 0,
        }
    print(f"[meta] Loaded {len(info)} virus accessions")
    return info


def load_virus_sequences(virus_dir, metadata):
    """加载病毒序列, 过滤掉 <1000 bp 和未分类的"""
    seqs = {}
    for f in sorted(Path(virus_dir).glob("*.fasta")):
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            if acc in metadata and metadata[acc]['length'] >= 1000:
                seqs[acc] = str(rec.seq)
    print(f"[seqs] Loaded {len(seqs)} sequences (>= 1000 bp)")
    return seqs


def generate_fragments(accessions, seqs, coverage_levels, n_per_cov, rng):
    """按覆盖率梯度为指定 accession 列表生成片段"""
    records = []
    idx = 0
    for cov_pct in coverage_levels:
        for _ in range(n_per_cov):
            for acc in accessions:
                if acc not in seqs: continue
                full_seq = seqs[acc]
                full_len = len(full_seq)
                frag_len = max(200, int(full_len * cov_pct / 100))
                frag_len = min(frag_len, full_len)
                max_start = full_len - frag_len
                start = rng.randint(0, max(1, max_start)) if max_start > 0 else 0
                frag = full_seq[start:start + frag_len]
                seq_id = f"family_leaveout|cov{cov_pct}|idx{idx:05d}|source={acc}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": full_len, "frag_length": len(frag),
                    "frag_seq": frag,
                })
                idx += 1
    return records


def generate_negatives(pos_lengths, host_genome, conserved_dir, eve_fasta,
                       n_A, n_B, n_C, rng):
    """生成三类负样本"""
    all_neg = []
    # A: host random
    host_seqs = []
    if host_genome and os.path.exists(host_genome):
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= 200: host_seqs.append(str(rec.seq))
    if not host_seqs: host_seqs = ["ACGT" * 5000]
    for i in range(n_A):
        target_len = rng.choice(pos_lengths) if pos_lengths else 500
        frag_len = max(200, target_len + rng.randint(-200, 200))
        seq = rng.choice(host_seqs)
        frag_len = min(frag_len, len(seq) - 5)
        start = rng.randint(0, max(1, len(seq) - frag_len))
        frag = seq[start:start + frag_len]
        all_neg.append({"seq_id": f"negative_A|negA_{i:04d}", "label": "negative_A",
                         "type": "negative_A", "frag_seq": frag, "frag_length": len(frag)})

    # B: conserved traps
    trap_seqs = []
    if conserved_dir and os.path.exists(conserved_dir):
        for f in Path(conserved_dir).glob("*.fasta"):
            try:
                rec = next(SeqIO.parse(f, "fasta"))
                if len(rec.seq) >= 500: trap_seqs.append(str(rec.seq))
            except: pass
    if not trap_seqs: trap_seqs = ["ACGT" * 3000]
    for i in range(n_B):
        seq = rng.choice(trap_seqs)
        frag_len = min(rng.choice(pos_lengths) if pos_lengths else 800, len(seq) - 5)
        frag_len = max(200, frag_len)
        start = rng.randint(0, max(1, len(seq) - frag_len))
        frag = seq[start:start + frag_len]
        all_neg.append({"seq_id": f"negative_B|negB_{i:04d}", "label": "negative_B",
                         "type": "negative_B", "frag_seq": frag, "frag_length": len(frag)})

    # C: EVE
    eve_seqs = []
    if eve_fasta and os.path.exists(eve_fasta):
        for rec in SeqIO.parse(eve_fasta, "fasta"):
            if len(rec.seq) >= 500: eve_seqs.append(str(rec.seq))
    if not eve_seqs:
        eve_seqs = host_seqs[:50] if host_seqs else ["ACGT" * 5000]
    for i in range(n_C):
        seq = rng.choice(eve_seqs)
        frag_len = min(rng.choice(pos_lengths) if pos_lengths else 1000, len(seq) - 5)
        frag_len = max(200, frag_len)
        start = rng.randint(0, max(1, len(seq) - frag_len))
        frag = seq[start:start + frag_len]
        all_neg.append({"seq_id": f"negative_C|negC_{i:04d}", "label": "negative_C",
                         "type": "negative_C", "frag_seq": frag, "frag_length": len(frag)})

    return all_neg


def write_dataset(outdir, pos_records, neg_records, prefix=""):
    """写出合并的FASTA和标签文件"""
    all_records = []
    all_labels = []
    for r in pos_records:
        all_records.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labels.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"],
                           "source": r.get("source", "")})
    for r in neg_records:
        all_records.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labels.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})

    fasta_out = os.path.join(outdir, f"{prefix}_eval.fasta")
    labels_out = os.path.join(outdir, f"{prefix}_labels.tsv")
    SeqIO.write(all_records, fasta_out, "fasta")
    pd.DataFrame(all_labels).to_csv(labels_out, sep='\t', index=False)
    return fasta_out, labels_out


def main():
    parser = argparse.ArgumentParser(description='按病毒科排除的未知病毒评估数据集')
    parser.add_argument('--virus-dir', required=True)
    parser.add_argument('--virus-meta', required=True, help='selected_viruses.tsv')
    parser.add_argument('--coverage-levels', type=int, nargs='+', default=[100, 80, 60, 40])
    parser.add_argument('--n-per-coverage', type=int, default=2)
    parser.add_argument('--host-genome')
    parser.add_argument('--conserved-prots')
    parser.add_argument('--eve-fasta')
    parser.add_argument('--n-neg-A', type=int, default=200)
    parser.add_argument('--n-neg-B', type=int, default=150)
    parser.add_argument('--n-neg-C', type=int, default=100)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    metadata = load_virus_metadata(args.virus_meta)
    seqs = load_virus_sequences(args.virus_dir, metadata)

    # 按科分组
    family_accessions = defaultdict(list)
    for acc, info in metadata.items():
        if acc in seqs:  # 仅 >=1000 bp
            family_accessions[info['family']].append(acc)

    # 排除只有一个病毒的科 (无法做留一法)
    families = {fam: accs for fam, accs in family_accessions.items()
                if len(accs) >= 2}
    families.pop('Unclassified', None)

    print(f"\n[family] {len(families)} families with ≥2 viruses:")
    for fam, accs in sorted(families.items(), key=lambda x: -len(x[1])):
        lens = [metadata[a]['length'] for a in accs]
        print(f"  {fam}: {len(accs)} viruses, {min(lens)}-{max(lens)} bp")

    # 预生成负样本 (所有科共用)
    pos_lengths_pool = [max(200, int(metadata[a]['length'] * c / 100))
                        for c in args.coverage_levels
                        for a in seqs]
    neg_records = generate_negatives(pos_lengths_pool, args.host_genome,
                                     args.conserved_prots, args.eve_fasta,
                                     args.n_neg_A, args.n_neg_B, args.n_neg_C, rng)
    # 只生成一份负样本 (所有科共用)
    neg_shared = neg_records[:args.n_neg_A + args.n_neg_B + args.n_neg_C]

    # 全局已知正样本 (所有科的病毒, 作为 "known" baseline)
    all_accs = list(seqs.keys())
    pos_all = generate_fragments(all_accs, seqs, args.coverage_levels, args.n_per_coverage, rng)
    write_dataset(args.outdir, pos_all, neg_shared, "all_known")

    # 每个科生成独立的 "留出" 数据集
    summary = []
    for fam, heldout_accs in sorted(families.items(), key=lambda x: -len(x[1])):
        fam_dir = os.path.join(args.outdir, f"leaveout_{fam.replace(' ', '_').replace('/','_')}")
        os.makedirs(fam_dir, exist_ok=True)

        pos_heldout = generate_fragments(heldout_accs, seqs, args.coverage_levels, args.n_per_coverage, rng)
        fasta, labels = write_dataset(fam_dir, pos_heldout, neg_shared, "heldout")

        # 其他科的病毒作为 "known" 正样本
        known_accs = [a for a in all_accs if a not in heldout_accs]
        pos_known = generate_fragments(known_accs, seqs, args.coverage_levels, args.n_per_coverage, rng)
        write_dataset(fam_dir, pos_known, neg_shared, "known")

        n_heldout = len(heldout_accs)
        summary.append({
            'family': fam,
            'n_viruses': n_heldout,
            'n_positive': len(pos_heldout),
            'mean_length': int(np.mean([metadata[a]['length'] for a in heldout_accs])),
        })
        print(f"  [{fam}] {n_heldout} viruses held out, {len(pos_heldout)} positive fragments")

    # 保存摘要
    pd.DataFrame(summary).to_csv(os.path.join(args.outdir, "family_leaveout_summary.tsv"), sep='\t', index=False)

    print(f"\n[DONE] {len(families)} family-leaveout datasets in {args.outdir}")
    print(f"  each with: heldout positives + known positives + shared negatives")
    print(f"  baseline: all_known_eval.fasta (all {len(all_accs)} viruses as known)")

if __name__ == '__main__':
    main()
