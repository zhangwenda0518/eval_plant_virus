#!/usr/bin/env python3
"""
构建候选病毒鉴定策略评估数据集（v4 — 简化版，无EVE）

正样本：从同一批50个评估用病毒按覆盖度梯度截取
  覆盖度: 100% 80% 60% 40% 20%，每个病毒每种覆盖度截取n条
  50病毒 × 5覆盖度 × n_per_cov = 250-500条

负样本A：从宿主基因组随机截取，长度分布与正样本匹配
负样本B：保守结构域陷阱序列

用法:
  python prep_build_id_eval_seqs.py \
      --virus-fasta eval_viruses_50/ \
      --coverage-levels 100 80 60 40 20 \
      --n-per-coverage 2 \
      --host-genome all.genome.uniq.fasta \
      --conserved-prots conserved_traps/ \
      --outdir eval_identification/ --seed 42
"""

import argparse, os, sys, random
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np


def generate_positive_by_coverage(virus_dir, out_dir, coverage_levels, n_per_cov, rng):
    """正样本：按覆盖度梯度从病毒基因组截取，随机起点模拟组装断点"""
    os.makedirs(out_dir, exist_ok=True)

    virus_seqs = []
    for f in sorted(Path(virus_dir).glob("*.fasta")):
        for rec in SeqIO.parse(f, "fasta"):
            if len(rec.seq) >= 200:
                virus_seqs.append((rec.id, str(rec.seq), len(rec.seq)))
    print(f"[positive] Loaded {len(virus_seqs)} virus sequences")

    if not virus_seqs:
        print("[ERROR] No virus sequences found!")
        return []

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

                seq_id = f"positive|cov{cov_pct}|idx{idx:05d}|source={acc}"
                rec = SeqRecord(Seq(frag), id=seq_id, description="")
                out_file = os.path.join(out_dir, f"pos_cov{cov_pct}_{idx:05d}.fasta")
                SeqIO.write(rec, out_file, "fasta")
                records.append({
                    "seq_id": seq_id, "label": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": full_len, "frag_length": frag_len,
                })
                idx += 1

    print(f"[positive] Generated {len(records)} sequences")
    return records


def load_host_sequences(host_genome, min_len=200):
    """加载宿主基因组"""
    seqs = []
    if host_genome and os.path.exists(host_genome):
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= min_len:
                seqs.append((rec.id, str(rec.seq)))
    return seqs


def generate_negative_A_length_matched(pos_records, host_seqs, out_dir, n, rng):
    """负样本A：从宿主基因组截取，长度与正样本匹配"""
    os.makedirs(out_dir, exist_ok=True)
    pos_lengths = [r["frag_length"] for r in pos_records]
    if not pos_lengths:
        pos_lengths = [rng.randint(200, 15000) for _ in range(n)]
    if not host_seqs:
        host_seqs = [("synthetic", "".join(rng.choices("ACGT", k=20000)))]

    records = []
    for i in range(n):
        target_len = rng.choice(pos_lengths)
        frag_len = max(200, target_len + rng.randint(-200, 200))

        chosen_id, chosen_seq = rng.choice(host_seqs)
        frag_len = min(frag_len, len(chosen_seq) - 5)
        if frag_len < 200:
            frag_len = min(200, len(chosen_seq))
        start = rng.randint(0, max(1, len(chosen_seq) - frag_len))
        frag = chosen_seq[start:start + frag_len]

        seq_id = f"negative_A|negA_{i:04d}"
        rec = SeqRecord(Seq(frag), id=seq_id, description="")
        out_file = os.path.join(out_dir, f"negA_{i:04d}.fasta")
        SeqIO.write(rec, out_file, "fasta")
        records.append({"seq_id": seq_id, "label": "negative_A", "frag_length": len(frag)})

    print(f"[negA] Generated {len(records)} sequences (length-matched)")
    return records


def generate_negative_B_traps(conserved_dir, out_dir, n, pos_records, rng):
    """负样本B：保守结构域陷阱序列，长度与正样本匹配"""
    os.makedirs(out_dir, exist_ok=True)
    pos_lengths = [r["frag_length"] for r in pos_records]

    # 尝试加载实际陷阱序列
    trap_seqs = []
    if conserved_dir and os.path.exists(conserved_dir):
        if os.path.isdir(conserved_dir):
            for f in Path(conserved_dir).glob("*.fasta"):
                recs = list(SeqIO.parse(f, "fasta"))
                if recs: trap_seqs.append((recs[0].id, str(recs[0].seq)))
        elif os.path.isfile(conserved_dir):
            for rec in SeqIO.parse(conserved_dir, "fasta"):
                if len(rec.seq) >= 500:
                    trap_seqs.append((rec.id, str(rec.seq)))

    records = []
    for i in range(n):
        if trap_seqs:
            chosen_id, chosen_seq = rng.choice(trap_seqs)
            frag_len = min(rng.choice(pos_lengths), len(chosen_seq) - 5)
            frag_len = max(200, frag_len)
            start = rng.randint(0, max(1, len(chosen_seq) - frag_len))
            frag = chosen_seq[start:start + frag_len]
            source = chosen_id
        else:
            # 无实际陷阱序列时生成合成序列（标记为placeholder）
            frag_len = rng.choice(pos_lengths)
            frag = "".join(rng.choices("ACGT", k=max(200, frag_len)))
            source = "synthetic"

        seq_id = f"negative_B|negB_{i:04d}|source={source}"
        rec = SeqRecord(Seq(frag), id=seq_id, description="")
        out_file = os.path.join(out_dir, f"negB_{i:04d}.fasta")
        SeqIO.write(rec, out_file, "fasta")
        records.append({"seq_id": seq_id, "label": "negative_B", "frag_length": len(frag)})

    print(f"[negB] Generated {len(records)} trap sequences ({len(trap_seqs)} real source)")
    return records


def load_eve_sequences(eve_fasta, host_genome):
    """
    加载 EVE/转座子序列。
    优先使用 --eve-fasta 提供的专门 EVE/转座子序列文件。
    降级方案：从宿主基因组中随机截取长序列（废弃了有Bug的氨基酸Motif扫描）。
    """
    eve_seqs = []
    if eve_fasta and os.path.exists(eve_fasta):
        for rec in SeqIO.parse(eve_fasta, "fasta"):
            if len(rec.seq) >= 500:
                eve_seqs.append((rec.id, str(rec.seq)))
        print(f"[eve] Loaded {len(eve_seqs)} EVE sequences from {eve_fasta}")
    elif host_genome and os.path.exists(host_genome):
        print("[eve] WARNING: No EVE FASTA provided. Falling back to sampling host genome for Neg-C.")
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= 5000:
                eve_seqs.append((rec.id, str(rec.seq)))
    return eve_seqs


def generate_negative_C_eve(eve_seqs, pos_records, out_dir, n, rng):
    """负样本C：从EVE/转座子序列截取，长度与正样本匹配。无序列时直接跳过（不生成合成序列）"""
    os.makedirs(out_dir, exist_ok=True)
    pos_lengths = [r["frag_length"] for r in pos_records]

    if not eve_seqs:
        print("[negC] ERROR: No sequences available to generate Negative C. Skipping.")
        return []

    records = []
    for i in range(n):
        chosen_id, chosen_seq = rng.choice(eve_seqs)
        target_len = rng.choice(pos_lengths) if pos_lengths else rng.randint(500, 5000)
        frag_len = min(target_len, len(chosen_seq) - 5)
        frag_len = max(200, frag_len)
        start = rng.randint(0, max(1, len(chosen_seq) - frag_len))
        frag = chosen_seq[start:start + frag_len)

        seq_id = f"negative_C|negC_{i:04d}|source={chosen_id}"
        rec = SeqRecord(Seq(frag), id=seq_id, description="")
        out_file = os.path.join(out_dir, f"negC_{i:04d}.fasta")
        SeqIO.write(rec, out_file, "fasta")
        records.append({"seq_id": seq_id, "label": "negative_C", "frag_length": len(frag)})

    print(f"[negC] Generated {len(records)} EVE/transposon sequences")
    return records


def merge_all(out_dirs, output_fasta, labels_tsv):
    """合并所有类别的序列为单一评估文件"""
    all_records = []
    all_labels = []

    for label, d in out_dirs.items():
        if not d or not os.path.exists(d):
            continue
        for f in Path(d).glob("*.fasta"):
            try:
                recs = list(SeqIO.parse(f, "fasta"))
                if not recs: continue
                rec = recs[0]
                all_records.append(rec)
                all_labels.append({"seq_id": rec.id, "label": label, "type": "positive" if "positive" in str(label) else "negative"})
            except Exception:
                continue

    SeqIO.write(all_records, output_fasta, "fasta")
    print(f"[merge] Wrote {len(all_records)} sequences to {output_fasta}")

    df = pd.DataFrame(all_labels)
    df.to_csv(labels_tsv, sep="\t", index=False)
    print(f"[merge] Wrote {len(df)} labels to {labels_tsv}")


def main():
    parser = argparse.ArgumentParser(description="构建鉴定评估数据集 (v4)")
    parser.add_argument("--virus-fasta", required=True, help="评估用病毒 FASTA 文件或目录")
    parser.add_argument("--coverage-levels", type=int, nargs="+",
                        default=[100, 80, 60, 40, 20], help="覆盖度梯度 (%)")
    parser.add_argument("--n-per-coverage", type=int, default=2,
                        help="每个病毒每种覆盖度的截取次数")
    parser.add_argument("--host-genome", help="宿主参考基因组 FASTA")
    parser.add_argument("--conserved-prots", help="保守结构域陷阱序列目录或FASTA")
    parser.add_argument("--eve-fasta", help="宿主EVE/转座子序列 FASTA (可选，无则从宿主基因组扫描RT模体)")
    parser.add_argument("--n-neg-A", type=int, default=500)
    parser.add_argument("--n-neg-B", type=int, default=300)
    parser.add_argument("--n-neg-C", type=int, default=200)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # 判断 virus-fasta 是目录还是文件
    virus_input = args.virus_fasta
    if os.path.isdir(virus_input):
        virus_dir = virus_input
    else:
        virus_dir = os.path.join(args.outdir, "_tmp_virus_split")
        os.makedirs(virus_dir, exist_ok=True)
        for rec in SeqIO.parse(virus_input, "fasta"):
            acc = rec.id.split()[0]
            SeqIO.write(rec, os.path.join(virus_dir, f"{acc}.fasta"), "fasta")

    out_dirs = {}

    # 正样本
    print("[1/3] Generating POSITIVE samples (coverage gradient)...")
    pos_dir = os.path.join(args.outdir, "positive")
    pos_records = generate_positive_by_coverage(virus_dir, pos_dir, args.coverage_levels, args.n_per_coverage, rng)
    out_dirs["positive"] = pos_dir

    # 宿主序列
    host_seqs = load_host_sequences(args.host_genome)

    # 负样本A
    print("[2/3] Generating NEGATIVE-A samples (host random, length-matched)...")
    negA_dir = os.path.join(args.outdir, "negative_A_random")
    negA_records = generate_negative_A_length_matched(pos_records, host_seqs, negA_dir, args.n_neg_A, rng)
    out_dirs["negative_A"] = negA_dir

    # 负样本B
    print("[3/4] Generating NEGATIVE-B samples (conserved domain traps)...")
    negB_dir = os.path.join(args.outdir, "negative_B_conserved")
    negB_records = generate_negative_B_traps(args.conserved_prots, negB_dir, args.n_neg_B, pos_records, rng)
    out_dirs["negative_B"] = negB_dir

    # 负样本C：EVE/转座子序列（宿主基因组中的内源性病毒元件）
    print("[4/4] Generating NEGATIVE-C samples (EVE/transposons)...")
    negC_dir = os.path.join(args.outdir, "negative_C_eve")
    eve_seqs = load_eve_sequences(args.eve_fasta, args.host_genome)
    negC_records = generate_negative_C_eve(eve_seqs, pos_records, negC_dir, args.n_neg_C, rng)
    out_dirs["negative_C"] = negC_dir

    # 合并输出
    merge_all(out_dirs,
              os.path.join(args.outdir, "evaluation_sequences.fasta"),
              os.path.join(args.outdir, "sequence_labels.tsv"))

    # 统计
    print(f"\n[SUMMARY]")
    print(f"  Positive: {len(pos_records)} ({len(args.coverage_levels)} coverages × {args.n_per_coverage} per × viruses)")
    print(f"  Negative A (host): {len(negA_records)}")
    print(f"  Negative B (traps): {len(negB_records)}")
    print(f"  Negative C (EVE/transposons): {len(negC_records)}")
    if pos_records:
        cov_dist = Counter(r["coverage_pct"] for r in pos_records)
        print(f"  Coverage distribution: {dict(sorted(cov_dist.items()))}")


if __name__ == "__main__":
    main()
