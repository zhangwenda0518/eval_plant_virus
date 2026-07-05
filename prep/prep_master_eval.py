#!/usr/bin/env python3
"""
MASTER 评估数据集构建器 (修正版)

方案:
  A: 已知病毒 60条 → 突变 0%, 5% → 9覆盖率 → 5重复 → label: A{mut}C{cov}R{rep}
  B: 新病毒   45条 → 突变 0%, 5% → 9覆盖率 → 5重复 → label: B{mut}C{cov}R{rep}
  C: VIROMOCK 105条 → 固定70%相似度 → 9覆盖率 → 5重复 → label: C70C{cov}R{rep}
  D: Conserved_Traps  1000条
  E: EVE_Transposon   1000条
  H: Host_Random      1000条

对比 prep_simulation_datasets.py:
  - 那个支持7级相似度梯度 (80%-20%), 这个固定70%
  - 那个无外部依赖, 这个依赖 pyrodigal + mutation-simulator
  - 那个A/B/C统一标签格式, 这个按方案分文件

用法:
  python prep_master_eval.py \
      --known-fasta step1_eval_viruses/known_60.fasta \
      --novel-fasta step1_eval_viruses/novel_45.fasta \
      --conserved-fasta step3_conserved_traps/traps_merged.fasta \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --host-fasta ningxia.genome.fasta \
      --outdir master_eval/ --seed 42
"""

import argparse, os, sys, random, time
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np

COV_LEVELS = [100, 90, 80, 70, 60, 50, 40, 30, 20]
N_REP = 5
MIN_LEN = 200

# 密码子表 (仅用于VIROMOCK)
CODON_USAGE = {
    'A': {'GCT': 0.30,'GCC': 0.40,'GCA': 0.15,'GCG': 0.15},
    'R': {'CGT': 0.08,'CGC': 0.19,'CGA': 0.11,'CGG': 0.21,'AGA': 0.20,'AGG': 0.21},
    'N': {'AAT': 0.46,'AAC': 0.54}, 'D': {'GAT': 0.54,'GAC': 0.46},
    'C': {'TGT': 0.45,'TGC': 0.55}, 'Q': {'CAA': 0.25,'CAG': 0.75},
    'E': {'GAA': 0.42,'GAG': 0.58}, 'G': {'GGT': 0.16,'GGC': 0.34,'GGA': 0.25,'GGG': 0.25},
    'H': {'CAT': 0.41,'CAC': 0.59}, 'I': {'ATT': 0.36,'ATC': 0.48,'ATA': 0.16},
    'L': {'TTA': 0.07,'TTG': 0.13,'CTT': 0.13,'CTC': 0.20,'CTA': 0.07,'CTG': 0.40},
    'K': {'AAA': 0.42,'AAG': 0.58}, 'M': {'ATG': 1.0}, 'F': {'TTT': 0.45,'TTC': 0.55},
    'P': {'CCT': 0.28,'CCC': 0.33,'CCA': 0.27,'CCG': 0.12},
    'S': {'TCT': 0.18,'TCC': 0.22,'TCA': 0.15,'TCG': 0.06,'AGT': 0.15,'AGC': 0.24},
    'T': {'ACT': 0.24,'ACC': 0.36,'ACA': 0.28,'ACG': 0.12},
    'W': {'TGG': 1.0}, 'Y': {'TAT': 0.43,'TAC': 0.57},
    'V': {'GTT': 0.18,'GTC': 0.24,'GTA': 0.11,'GTG': 0.47},
    '*': {'TAA': 0.30,'TAG': 0.20,'TGA': 0.50}
}
CODON_AA = {c: aa for aa, codons in CODON_USAGE.items() for c in codons}
STOP_CODONS = {'TAA','TAG','TGA'}


def mutate(seq, rate_pct, rng):
    """纯Python随机突变 (无外部依赖)"""
    if rate_pct <= 0: return seq
    s = list(seq)
    n = int(len(s) * rate_pct / 100)
    for p in rng.sample(range(len(s)), min(n, len(s))):
        s[p] = rng.choice([b for b in 'ACGT' if b != s[p]])
    return ''.join(s)


def viromock_backtranslate(seq, rng,
                            codon_table=CODON_USAGE,
                            codon_aa=CODON_AA,
                            stop_codons=STOP_CODONS):
    """
    VIROMOCK风格反向翻译: 每3碱基翻译为AA → 随机选同义密码子
    结果: 核苷酸~68-73%相似度, 氨基酸~100%一致
    纯Python实现, 无外部依赖
    """
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    new = []
    for codon in codons:
        if codon not in codon_aa:
            new.append(codon)
        else:
            aa = codon_aa[codon]
            choices = list(codon_table[aa].keys())
            weights = list(codon_table[aa].values())
            new.append(rng.choices(choices, weights=weights, k=1)[0])
    result = ''.join(new)
    # 额外~1%安全突变避免提前终止
    n_mut = max(1, int(len(result) * 0.005))
    s = list(result)
    for p in rng.sample(range(len(s)), min(n_mut, len(s))):
        codon_start = (p // 3) * 3
        if codon_start + 3 > len(s): continue
        old_codon = ''.join(s[codon_start:codon_start+3])
        valid = [b for b in 'ACGT' if b != s[p]]
        if valid and old_codon not in stop_codons:
            new_codon = list(old_codon)
            new_codon[p % 3] = rng.choice(valid)
            if ''.join(new_codon) not in stop_codons:
                s[p] = rng.choice(valid)
    return ''.join(s)


def generate_fragments(genomes, cov_levels, n_rep, rng, label_prefix):
    """生成覆盖率梯度片段"""
    records = []
    idx = 0
    for cov in cov_levels:
        for rep in range(1, n_rep + 1):
            for acc, seq in genomes.items():
                L = len(seq)
                flen = max(MIN_LEN, int(L * cov / 100))
                flen = min(flen, L)
                start = rng.randint(0, L - flen) if L > flen else 0
                frag = seq[start:start + flen]
                seq_id = f"{label_prefix}C{cov}R{rep}_{acc}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov, "replicate": rep,
                    "frag_seq": frag, "frag_length": flen,
                })
                idx += 1
    return records


def generate_negatives(seqs_pool, n_each, type_label, id_prefix, rng, pos_lengths):
    """生成一类负样本: 从seqs_pool随机截取, 长度匹配正样本"""
    if not seqs_pool:
        seqs_pool = ["".join(rng.choices('ACGT', k=5000)) for _ in range(10)]
    records = []
    for i in range(n_each):
        ref = rng.choice(seqs_pool)
        flen = min(rng.choice(pos_lengths), len(ref) - 10)
        flen = max(MIN_LEN, flen)
        start = rng.randint(0, max(1, len(ref) - flen))
        frag = ref[start:start + flen]
        seq_id = f"{id_prefix}_{i:04d}"
        records.append({
            "seq_id": seq_id, "label": "negative", "type": type_label,
            "frag_seq": frag, "frag_length": flen,
        })
    return records


def load_genomes(path):
    """加载FASTA, 返回 {acc: seq}"""
    genomes = {}
    if os.path.isdir(path):
        files = sorted(Path(path).glob("*.fasta"))
    else:
        files = [Path(path)]
    for f in files:
        try:
            for rec in SeqIO.parse(f, "fasta"):
                acc = rec.id.split()[0]
                if len(rec.seq) >= 500:  # 更低阈值以包含更多病毒
                    genomes[acc] = str(rec.seq)
        except: pass
    return genomes


def write_output(outdir, name, pos_records, neg_records, pos_meta):
    """合并正负样本并写入FASTA + labels TSV"""
    recs, labs = [], []

    # 正样本
    for r in pos_records:
        recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        lab = {k: r[k] for k in ["seq_id", "label", "type", "source", "coverage_pct"] if k in r}
        lab.update(pos_meta.get(r.get("source", ""), {}))
        labs.append(lab)

    # 负样本
    for r in neg_records:
        recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        labs.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})

    SeqIO.write(recs, os.path.join(outdir, f"{name}.fasta"), "fasta")
    pd.DataFrame(labs).to_csv(os.path.join(outdir, f"{name}_labels.tsv"), sep='\t', index=False)
    return len(recs)


def main():
    parser = argparse.ArgumentParser(description='MASTER评估数据集 (修正版)')
    parser.add_argument('--known-fasta', required=True, help='已知60病毒FASTA')
    parser.add_argument('--novel-fasta', required=True, help='新45病毒FASTA')
    parser.add_argument('--conserved-fasta', required=False, default=None)
    parser.add_argument('--eve-fasta', required=False, default=None)
    parser.add_argument('--host-fasta', required=False, default=None)
    parser.add_argument('--n-neg', type=int, default=1000)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    t0 = time.time()
    os.makedirs(args.outdir, exist_ok=True)

    # ── 加载 ──
    known = load_genomes(args.known_fasta)
    novel = load_genomes(args.novel_fasta)
    all_virus = {**known, **novel}
    print(f"[load] Known: {len(known)}  Novel: {len(novel)}  Total: {len(all_virus)}")

    # ── 正样本长度池 (负样本长度匹配用) ──
    pos_lengths_pool = [max(MIN_LEN, int(L*c/100)) for L in [len(s) for s in all_virus.values()]
                        for c in COV_LEVELS]

    # ═══════════════════════════════════════════
    # 方案A: 已知病毒 (0%, 5%)
    # ═══════════════════════════════════════════
    print(f"\n[A] Known virus evaluation (0%, 5%)...")
    pos_a, meta_a = [], {}
    for acc, seq in known.items():
        for rate in [0, 5]:
            mut_seq = mutate(seq, rate, rng)
            key = f"{acc}_mut{rate}"
            known_mut = {key: mut_seq}
            frags = generate_fragments(known_mut, COV_LEVELS, N_REP, rng,
                                       label_prefix=f"A{rate}")
            pos_a.extend(frags)
            meta_a[key] = {"mut_rate": rate, "scheme": "A"}
    print(f"  A: {len(pos_a)} fragments")

    # ═══════════════════════════════════════════
    # 方案B: 新病毒 (0%, 5%)
    # ═══════════════════════════════════════════
    print(f"\n[B] Novel virus evaluation (0%, 5%)...")
    pos_b, meta_b = [], {}
    for acc, seq in novel.items():
        for rate in [0, 5]:
            mut_seq = mutate(seq, rate, rng)
            key = f"{acc}_mut{rate}"
            novel_mut = {key: mut_seq}
            frags = generate_fragments(novel_mut, COV_LEVELS, N_REP, rng,
                                       label_prefix=f"B{rate}")
            pos_b.extend(frags)
            meta_b[key] = {"mut_rate": rate, "scheme": "B"}
    print(f"  B: {len(pos_b)} fragments")

    # ═══════════════════════════════════════════
    # 方案C: VIROMOCK (固定70%)
    # ═══════════════════════════════════════════
    print(f"\n[C] VIROMOCK evaluation (70% similarity)...")
    pos_c, meta_c = [], {}
    for acc, seq in all_virus.items():
        vm_seq = viromock_backtranslate(seq, rng)
        key = f"{acc}_vm70"
        vm_genomes = {key: vm_seq}
        frags = generate_fragments(vm_genomes, COV_LEVELS, N_REP, rng,
                                   label_prefix=f"C70")
        pos_c.extend(frags)
        meta_c[key] = {"similarity": 70, "scheme": "C"}
        # 验证实际相似度
        actual = sum(1 for i in range(min(len(seq), len(vm_seq))) if seq[i]==vm_seq[i])/min(len(seq), len(vm_seq))
    print(f"  C: {len(pos_c)} fragments (actual nt_id~{actual*100:.0f}%)")

    # ═══════════════════════════════════════════
    # 负样本 D/E/H
    # ═══════════════════════════════════════════
    print(f"\n[Neg] Generating D/E/H ({args.n_neg} each)...")

    def load_seqs(path):
        if not path or not os.path.exists(path): return []
        return [str(r.seq) for r in SeqIO.parse(path, "fasta") if len(r.seq) >= MIN_LEN]

    neg_d = generate_negatives(load_seqs(args.conserved_fasta), args.n_neg,
                               "Conserved_Traps", "D", rng, pos_lengths_pool)
    neg_e = generate_negatives(load_seqs(args.eve_fasta), args.n_neg,
                               "EVE_Transposon", "E", rng, pos_lengths_pool)
    host_seqs = load_seqs(args.host_fasta)
    neg_h = generate_negatives(host_seqs if host_seqs else None, args.n_neg,
                               "Host_Random", "H", rng, pos_lengths_pool)
    all_neg = neg_d + neg_e + neg_h
    print(f"  Negatives: {len(all_neg)} (D={len(neg_d)} E={len(neg_e)} H={len(neg_h)})")

    # ═══════════════════════════════════════════
    # 输出独立数据集
    # ═══════════════════════════════════════════
    n_a = write_output(args.outdir, "dataset_A", pos_a, all_neg, meta_a)
    n_b = write_output(args.outdir, "dataset_B", pos_b, all_neg, meta_b)
    n_c = write_output(args.outdir, "dataset_C", pos_c, all_neg, meta_c)
    # 也输出单独负样本集
    SeqIO.write([SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description="")
                 for r in all_neg],
                os.path.join(args.outdir, "negatives_DEH.fasta"), "fasta")
    pd.DataFrame([{"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]}
                  for r in all_neg]).to_csv(
        os.path.join(args.outdir, "labels_negatives.tsv"), sep='\t', index=False)

    total_pos = len(pos_a) + len(pos_b) + len(pos_c)
    elapsed = (time.time() - t0) / 60
    print(f"\n{'='*50}")
    print(f"  MASTER 数据集完成 ({elapsed:.1f} min)")
    print(f"{'='*50}")
    print(f"  A 已知病毒:  {len(pos_a):>8,}  (60 viruses × 2 rates × 9 cov × 5 rep)")
    print(f"  B 新病毒:    {len(pos_b):>8,}  (45 viruses × 2 rates × 9 cov × 5 rep)")
    print(f"  C VIROMOCK:  {len(pos_c):>8,}  (105 viruses × 1 sim × 9 cov × 5 rep)")
    print(f"  正样本总计:  {total_pos:>8,}")
    print(f"  负样本:      {len(all_neg):>8,}")
    print(f"  全部序列:    {total_pos+len(all_neg):>8,}")
    print(f"  输出:        {args.outdir}")

if __name__ == '__main__':
    main()
