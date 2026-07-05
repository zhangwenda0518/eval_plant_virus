#!/usr/bin/env python3
"""
统一的新病毒鉴定评估数据集构建器
策略: 四层未知性金字塔 (留一科 + 突变梯度)

层0 (已知基线):  0% 突变,  参考库有该病毒
层1 (株系变异):  5-10% 突变, 模拟同种不同分离物
层2 (近缘新种):  15-25% 突变, 模拟同属未收录物种
层3 (远缘新属):  30-40% 突变, 模拟同科未收录属
层4 (全新病毒科): 留一科, 模拟完全没有近缘参考的全新病毒

用法:
  python prep_build_novel_virus_unified.py \
      --virus-dir step1_eval_viruses/ \
      --virus-meta step1_eval_viruses/selected_viruses.tsv \
      --coverage-levels 100 80 60 40 --n-per-coverage 2 \
      --host-genome ningxia.genome.fasta \
      --conserved-prots step3_conserved_traps/ \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --outdir step3_novel_unified/ --seed 42
"""

import argparse, os, sys, random, subprocess, shutil
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np


# ── 数据加载 ──

def load_virus_metadata(meta_path):
    meta = pd.read_csv(meta_path, sep='\t')
    info = {}
    for _, row in meta.iterrows():
        acc = row['accession']
        info[acc] = {
            'family': str(row['family']), 'genus': str(row['genus']),
            'species': str(row['species']), 'genome_type': str(row.get('genome_type', '')),
            'length': int(row['length'])
        }
    return info

def load_virus_sequences(virus_dir, metadata, min_len=1000):
    seqs = {}
    for f in Path(virus_dir).glob("*.fasta"):
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            if acc in metadata and metadata[acc]['length'] >= min_len:
                seqs[acc] = str(rec.seq)
    return seqs


# ── 突变生成 ──

def mutate_sequence(seq, rate, rng):
    """纯Python随机突变, 无需外部工具"""
    if rate == 0: return seq
    seq_list = list(seq)
    n_mut = int(len(seq) * rate / 100)
    positions = rng.sample(range(len(seq)), n_mut)
    bases = ['A', 'C', 'G', 'T']
    for pos in positions:
        original = seq_list[pos]
        choices = [b for b in bases if b != original]
        seq_list[pos] = rng.choice(choices) if choices else original
    return ''.join(seq_list)


# ── 片段生成 ──

def generate_fragments(accessions, seqs, coverage_levels, n_per_cov, rng,
                       label_prefix="novel", extra_meta=None):
    """按覆盖率梯度截取片段"""
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

                # 序列ID编码所有元信息
                meta_str = ""
                if extra_meta:
                    meta_str = '|' + '|'.join(f'{k}={v}' for k, v in extra_meta.items())
                seq_id = f"{label_prefix}|cov{cov_pct}|idx{idx:05d}|source={acc}{meta_str}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": full_len, "frag_length": len(frag),
                    "frag_seq": frag,
                })
                idx += 1
    return records


# ── 负样本 ──

def generate_negatives(pos_lengths, host_genome, conserved_dir, eve_fasta,
                       n_A, n_B, n_C, rng):
    """生成三类负样本"""
    # A: host random
    host_seqs = []
    if host_genome and os.path.exists(host_genome):
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= 200: host_seqs.append(str(rec.seq))
    if not host_seqs: host_seqs = ["ACGT" * 5000]

    neg_A = []
    for i in range(n_A):
        tlen = rng.choice(pos_lengths) if pos_lengths else 500
        flen = max(200, tlen + rng.randint(-200, 200))
        seq = rng.choice(host_seqs)
        flen = min(flen, len(seq) - 5)
        start = rng.randint(0, max(1, len(seq) - flen))
        neg_A.append({"seq_id": f"negA|{i:04d}", "label": "negative", "type": "negative_A",
                       "frag_seq": seq[start:start+flen], "frag_length": flen})

    # B: conserved traps
    trap_seqs = []
    if conserved_dir and os.path.exists(conserved_dir):
        for f in Path(conserved_dir).glob("*.fasta"):
            try:
                rec = next(SeqIO.parse(f, "fasta"))
                if len(rec.seq) >= 500: trap_seqs.append(str(rec.seq))
            except: pass
    if not trap_seqs: trap_seqs = ["ACGT" * 3000]
    neg_B = []
    for i in range(n_B):
        seq = rng.choice(trap_seqs)
        flen = min(rng.choice(pos_lengths) if pos_lengths else 800, len(seq) - 5)
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(seq) - flen))
        neg_B.append({"seq_id": f"negB|{i:04d}", "label": "negative", "type": "negative_B",
                       "frag_seq": seq[start:start+flen], "frag_length": flen})

    # C: EVE
    eve_seqs = []
    if eve_fasta and os.path.exists(eve_fasta):
        for rec in SeqIO.parse(eve_fasta, "fasta"):
            if len(rec.seq) >= 500: eve_seqs.append(str(rec.seq))
    if not eve_seqs: eve_seqs = host_seqs[:50]
    neg_C = []
    for i in range(n_C):
        seq = rng.choice(eve_seqs)
        flen = min(rng.choice(pos_lengths) if pos_lengths else 1000, len(seq) - 5)
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(seq) - flen))
        neg_C.append({"seq_id": f"negC|{i:04d}", "label": "negative", "type": "negative_C",
                       "frag_seq": seq[start:start+flen], "frag_length": flen})

    return neg_A + neg_B + neg_C


def write_dataset(outdir, name, pos_records, neg_records):
    fasta_out = os.path.join(outdir, f"{name}.fasta")
    labels_out = os.path.join(outdir, f"{name}.tsv")
    all_recs, all_labs = [], []
    for r in pos_records:
        all_recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labs.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"],
                         "source": r.get("source", ""), "coverage_pct": r.get("coverage_pct","")})
    for r in neg_records:
        all_recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labs.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})
    SeqIO.write(all_recs, fasta_out, "fasta")
    pd.DataFrame(all_labs).to_csv(labels_out, sep='\t', index=False)
    return fasta_out


# ── 主函数 ──

def main():
    parser = argparse.ArgumentParser(description='统一四层新病毒评估数据集')
    parser.add_argument('--virus-dir', required=True)
    parser.add_argument('--virus-meta', required=True)
    parser.add_argument('--coverage-levels', type=int, nargs='+', default=[100, 80, 60, 40])
    parser.add_argument('--n-per-coverage', type=int, default=2)
    parser.add_argument('--host-genome')
    parser.add_argument('--conserved-prots')
    parser.add_argument('--eve-fasta')
    parser.add_argument('--n-neg', type=int, default=300, help='每类负样本数量')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    metadata = load_virus_metadata(args.virus_meta)
    seqs = load_virus_sequences(args.virus_dir, metadata, min_len=1000)
    all_accs = list(seqs.keys())
    print(f"[init] {len(all_accs)} viruses >= 1000 bp")

    # 按科分组
    fam2accs = defaultdict(list)
    for acc in all_accs:
        fam2accs[metadata[acc]['family']].append(acc)
    families = {f: a for f, a in fam2accs.items() if len(a) >= 2 and f != 'Unclassified'}
    # 选前5大科做突变梯度, 其余只用留出法
    top_families = sorted(families.items(), key=lambda x: -len(x[1]))[:5]
    top_fam_names = {f for f, _ in top_families}
    print(f"[family] {len(families)} families with ≥2 viruses")
    print(f"[family] Top 5 for mutation gradient: {[f for f,_ in top_families]}")

    # 预生成共享负样本
    pos_lens = [max(200, int(metadata[a]['length'] * c / 100)) for a in all_accs
                for c in args.coverage_levels]
    neg_records = generate_negatives(pos_lens, args.host_genome, args.conserved_prots,
                                     args.eve_fasta, args.n_neg, args.n_neg, args.n_neg, rng)

    # ═══════════════════════════════════════════
    # 层0: 已知基线 (所有病毒, 0%突变)
    # ═══════════════════════════════════════════
    print("\n[Layer 0] Known baseline (0% mutation)...")
    pos_l0 = generate_fragments(all_accs, seqs, args.coverage_levels, args.n_per_coverage, rng,
                                label_prefix="known", extra_meta={"mut": 0, "layer": 0})
    write_dataset(args.outdir, "layer0_known", pos_l0, neg_records)

    # ═══════════════════════════════════════════
    # 层1-3: 突变梯度 (针对5个最大科)
    # ═══════════════════════════════════════════
    mut_layers = {
        1: {'rates': [5, 10], 'label': 'strain_variant'},
        2: {'rates': [15, 25], 'label': 'novel_species'},
        3: {'rates': [30, 40], 'label': 'novel_genus'},
    }

    all_mut_pos = []
    for layer, cfg in mut_layers.items():
        for rate in cfg['rates']:
            label = cfg['label']
            print(f"\n[Layer {layer}] {label} ({rate}% mutation)...")
            mut_seqs = {}
            for acc in all_accs:
                fam = metadata[acc]['family']
                if fam in top_fam_names:
                    mut_seqs[f"{acc}_mut{rate}pct"] = mutate_sequence(seqs[acc], rate, rng)

            pos = generate_fragments(list(mut_seqs.keys()), mut_seqs,
                                     args.coverage_levels, args.n_per_coverage, rng,
                                     label_prefix=label,
                                     extra_meta={"mut": rate, "layer": layer})
            all_mut_pos.extend(pos)

    write_dataset(args.outdir, "layer1-3_mutation_gradient", all_mut_pos, neg_records)

    # ═══════════════════════════════════════════
    # 层4: 留一科交叉验证 (全14个科)
    # ═══════════════════════════════════════════
    print(f"\n[Layer 4] Leave-one-family-out ({len(families)} families)...")
    summary = []
    for fam, heldout_accs in sorted(families.items(), key=lambda x: -len(x[1])):
        fam_safe = fam.replace(' ', '_').replace('/', '_')
        fam_dir = os.path.join(args.outdir, f"layer4_leaveout_{fam_safe}")
        os.makedirs(fam_dir, exist_ok=True)

        # 该科的病毒序列
        pos = generate_fragments(heldout_accs, seqs, args.coverage_levels, args.n_per_coverage, rng,
                                 label_prefix="novel_family",
                                 extra_meta={"heldout_family": fam, "layer": 4})
        write_dataset(fam_dir, f"heldout_{fam_safe}", pos, neg_records)

        # 其他科作为 "known" baseline
        known_accs = [a for a in all_accs if a not in heldout_accs]
        pos_known = generate_fragments(known_accs, seqs, args.coverage_levels, args.n_per_coverage, rng,
                                       label_prefix="known", extra_meta={"layer": 0})
        write_dataset(fam_dir, f"known_context", pos_known, neg_records)

        n_heldout_seqs = len(pos)
        summary.append({
            'family': fam, 'n_heldout_viruses': len(heldout_accs),
            'n_heldout_fragments': n_heldout_seqs,
            'mean_length': int(np.mean([metadata[a]['length'] for a in heldout_accs])),
            'genome_types': ','.join(set(metadata[a]['genome_type'] for a in heldout_accs)),
        })
        print(f"  [{fam}] {len(heldout_accs)} viruses → {n_heldout_seqs} fragments")

    pd.DataFrame(summary).to_csv(os.path.join(args.outdir, "layer4_family_summary.tsv"), sep='\t', index=False)

    # ═══════════════════════════════════════════
    # 汇总清单
    # ═══════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"  四层新病毒评估数据集生成完毕")
    print(f"{'='*60}")
    print(f"  层0 (已知基线):           layer0_known.fasta")
    print(f"  层1-3 (突变梯度):         layer1-3_mutation_gradient.fasta")
    print(f"  层4 (留一科):             layer4_leaveout_*/ (共{len(families)}个科)")
    print(f"  共享负样本:               {args.n_neg*3} 条 (A/B/C 各{args.n_neg})")
    print(f"  输出目录:                 {args.outdir}")

if __name__ == '__main__':
    main()
