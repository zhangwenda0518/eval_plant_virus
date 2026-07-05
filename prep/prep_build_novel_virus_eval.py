#!/usr/bin/env python3
"""
构建未知病毒鉴定评估数据集 (Novel Virus Identification Evaluation)

策略:
  1. 选取 N 个代表性病毒作为"未知病毒"候选
  2. 用 virome_simulator.py mutate 生成 5 个突变梯度 (0%–40%)
  3. 按覆盖率梯度截取片段 (同 prep_build_id_eval_seqs.py)
  4. 加入负样本 (宿主随机 + 保守结构域陷阱 + EVE)
  5. 输出合并评估文件

用法:
  python prep_build_novel_virus_eval.py \
      --virus-dir step1_eval_viruses/ \
      --n-viruses 10 \
      --mut-rates 0 10 20 30 40 \
      --coverage-levels 100 80 60 40 \
      --n-per-coverage 2 \
      --host-genome ningxia.genome.fasta \
      --conserved-prots step3_conserved_traps/ \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --n-neg-A 200 --n-neg-B 150 --n-neg-C 100 \
      --outdir step3_novel_virus_eval/ --seed 42
"""

import argparse, os, sys, random, subprocess, shutil, tempfile
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np


def select_novel_viruses(virus_dir, n, rng):
    """从评估病毒中选取 n 个作为未知病毒 (优先选不同科属)"""
    seqs = []
    for f in sorted(Path(virus_dir).glob("*.fasta")):
        for rec in SeqIO.parse(f, "fasta"):
            if len(rec.seq) >= 1000:  # 排除类病毒
                seqs.append((rec.id, str(rec.seq), len(rec.seq)))

    # 简单随机选取 (如果有 metadata 可按科属分层)
    selected = rng.sample(seqs, min(n, len(seqs)))
    print(f"[select] {len(selected)} viruses for novel evaluation:")
    for acc, _, length in selected:
        print(f"  {acc}: {length} bp")
    return selected


def mutate_genomes(virus_seqs, rates, outdir, rng):
    """调用 mutation-simulator 生成突变基因组 (或直接用 Python 随机突变)"""
    os.makedirs(outdir, exist_ok=True)
    mutated = {}  # rate -> [(acc, seq, length), ...]

    for rate in rates:
        rate_dir = os.path.join(outdir, f"mut_{rate}pct")
        os.makedirs(rate_dir, exist_ok=True)
        mut_seqs = []

        for acc, seq, length in virus_seqs:
            if rate == 0:
                mut_seq = seq
            else:
                # 简单随机突变 (无需外部工具)
                mut_seq = list(seq)
                n_mut = int(length * rate / 100)
                positions = rng.sample(range(length), n_mut)
                bases = ['A', 'C', 'G', 'T']
                for pos in positions:
                    original = mut_seq[pos]
                    mut_seq[pos] = rng.choice([b for b in bases if b != original])
                mut_seq = ''.join(mut_seq)

            out_path = os.path.join(rate_dir, f"{acc.split('.')[0]}_mut{rate}pct.fasta")
            rec = SeqRecord(Seq(mut_seq), id=f"{acc}_mut{rate}pct", description=f"mutation_rate={rate}%")
            SeqIO.write(rec, out_path, "fasta")
            mut_seqs.append((rec.id, mut_seq, length))

        mutated[rate] = mut_seqs
        print(f"[mutate] {rate}%: {len(mut_seqs)} genomes in {rate_dir}")

    return mutated


def generate_fragments(virus_seqs, coverage_levels, n_per_cov, rng):
    """按覆盖率梯度截取片段 (同 prep_build_id_eval_seqs.py)"""
    records = []
    idx = 0
    for cov_pct in coverage_levels:
        for _ in range(n_per_cov):
            for acc, full_seq, full_len in virus_seqs:
                frag_len = max(200, int(full_len * cov_pct / 100))
                frag_len = min(frag_len, full_len)
                max_start = full_len - frag_len
                start = rng.randint(0, max(1, max_start)) if max_start > 0 else 0
                frag = full_seq[start:start + frag_len]
                seq_id = f"novel|mut{acc.split('_mut')[-1].replace('pct','')}pct|cov{cov_pct}|idx{idx:05d}|source={acc}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": full_len, "frag_length": len(frag),
                    "frag_seq": frag,
                })
                idx += 1
    return records


def generate_negatives(pos_records, host_genome, conserved_dir, eve_fasta,
                       n_neg_A, n_neg_B, n_neg_C, rng):
    """生成三类负样本 (简化版, 复用 prep_build_id_eval_seqs 逻辑)"""
    all_neg = []
    pos_lengths = [r["frag_length"] for r in pos_records]

    # Negative A: host random
    host_seqs = []
    if host_genome and os.path.exists(host_genome):
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= 200: host_seqs.append((rec.id, str(rec.seq)))
    if not host_seqs: host_seqs = [("syn", "ACGT" * 5000)]
    for i in range(n_neg_A):
        target_len = rng.choice(pos_lengths)
        frag_len = max(200, target_len + rng.randint(-200, 200))
        _, chosen_seq = rng.choice(host_seqs)
        frag_len = min(frag_len, len(chosen_seq) - 5)
        start = rng.randint(0, max(1, len(chosen_seq) - frag_len))
        frag = chosen_seq[start:start + frag_len]
        all_neg.append({"seq_id": f"negative_A|negA_{i:04d}", "label": "negative_A",
                         "type": "negative_A", "frag_seq": frag, "frag_length": len(frag)})

    # Negative B: conserved traps
    trap_seqs = []
    if conserved_dir and os.path.exists(conserved_dir):
        for f in Path(conserved_dir).glob("*.fasta"):
            try:
                rec = next(SeqIO.parse(f, "fasta"))
                if len(rec.seq) >= 500: trap_seqs.append((rec.id, str(rec.seq)))
            except: pass
    if not trap_seqs: trap_seqs = [("syn", "ACGT" * 3000)]
    for i in range(n_neg_B):
        _, chosen_seq = rng.choice(trap_seqs)
        frag_len = min(rng.choice(pos_lengths), len(chosen_seq) - 5)
        frag_len = max(200, frag_len)
        start = rng.randint(0, max(1, len(chosen_seq) - frag_len))
        frag = chosen_seq[start:start + frag_len]
        all_neg.append({"seq_id": f"negative_B|negB_{i:04d}", "label": "negative_B",
                         "type": "negative_B", "frag_seq": frag, "frag_length": len(frag)})

    # Negative C: EVE
    eve_seqs = []
    if eve_fasta and os.path.exists(eve_fasta):
        for rec in SeqIO.parse(eve_fasta, "fasta"):
            if len(rec.seq) >= 500: eve_seqs.append((rec.id, str(rec.seq)))
    if not eve_seqs:
        # fallback to host genome
        if host_seqs:
            eve_seqs = [(f"host_{i}", s) for i, (_, s) in enumerate(host_seqs[:50]) if len(s) >= 5000]
    if eve_seqs:
        for i in range(n_neg_C):
            _, chosen_seq = rng.choice(eve_seqs)
            frag_len = min(rng.choice(pos_lengths), len(chosen_seq) - 5)
            frag_len = max(200, frag_len)
            start = rng.randint(0, max(1, len(chosen_seq) - frag_len))
            frag = chosen_seq[start:start + frag_len]
            all_neg.append({"seq_id": f"negative_C|negC_{i:04d}", "label": "negative_C",
                             "type": "negative_C", "frag_seq": frag, "frag_length": len(frag)})

    return all_neg


def main():
    parser = argparse.ArgumentParser(description='构建未知病毒鉴定评估数据集')
    parser.add_argument('--virus-dir', required=True, help='step1_eval_viruses/')
    parser.add_argument('--n-viruses', type=int, default=10, help='选取的未知病毒数量')
    parser.add_argument('--mut-rates', type=int, nargs='+', default=[0, 10, 20, 30, 40])
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

    # 1. 选取未知病毒
    selected = select_novel_viruses(args.virus_dir, args.n_viruses, rng)

    # 2. 突变
    mutated = mutate_genomes(selected, args.mut_rates, args.outdir, rng)

    # 3. 生成正样本 (每个突变率 × 覆盖率 × 复制)
    all_pos = []
    for rate in args.mut_rates:
        pos = generate_fragments(mutated[rate], args.coverage_levels, args.n_per_coverage, rng)
        all_pos.extend(pos)
    print(f"[positive] Total: {len(all_pos)} fragments "
          f"({args.n_viruses} viruses × {len(args.mut_rates)} rates × {len(args.coverage_levels)} cov × {args.n_per_coverage} rep)")

    # 4. 生成负样本
    all_neg = generate_negatives(all_pos, args.host_genome, args.conserved_prots,
                                 args.eve_fasta, args.n_neg_A, args.n_neg_B, args.n_neg_C, rng)

    # 5. 合并输出
    all_records = []
    all_labels = []
    for r in all_pos:
        all_records.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labels.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})
    for r in all_neg:
        all_records.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labels.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})

    fasta_out = os.path.join(args.outdir, "novel_virus_eval.fasta")
    labels_out = os.path.join(args.outdir, "novel_virus_labels.tsv")
    SeqIO.write(all_records, fasta_out, "fasta")
    pd.DataFrame(all_labels).to_csv(labels_out, sep='\t', index=False)

    # 统计
    pos_count = len(all_pos)
    neg_count = len(all_neg)
    print(f"\n[DONE]")
    print(f"  Total sequences: {pos_count + neg_count}")
    print(f"  Positive: {pos_count} (novel virus, 5 mutation rates)")
    print(f"  Negative: {neg_count} (A={args.n_neg_A} B={args.n_neg_B} C={args.n_neg_C})")
    print(f"  FASTA: {fasta_out}")
    print(f"  Labels: {labels_out}")

if __name__ == '__main__':
    main()
