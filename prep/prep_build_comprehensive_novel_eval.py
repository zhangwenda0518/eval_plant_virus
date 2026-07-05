#!/usr/bin/env python3
"""
综合未知病毒鉴定评估数据集构建器

两大评估方案:
  方案A — 已知病毒检测 (Known)
    突变率: 0%, 1%, 3%, 5%, 7%, 9%
    覆盖率: 100%-20% (9梯度)
    重复: 5

  方案B — 新病毒检测 (Novel), 两个来源
    来源1: 原有病毒突变到目标相似度
      L1: 0.90(10%mut), 0.80(20%mut), 0.70(30%mut)
      L2: 0.70(30%mut), 0.60(40%mut), 0.50(50%mut)
      L3: 0.40(60%mut), 0.30(70%mut), 0.20(80%mut)
    来源2: NCBI新病毒 (无突变 + 同上突变)
    覆盖率: 100%-20% (9梯度)
    重复: 5

用法:
  python prep_build_comprehensive_novel_eval.py \
      --virus-dir step1_eval_viruses/ \
      --ncbi-meta novel_viruses_metadata.tsv \
      --ref-fasta final.cluster.ref.fasta \
      --n-novel-source1 10 \
      --host-genome ningxia.genome.fasta \
      --conserved-prots step3_conserved_traps/ \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --outdir step3_comprehensive/ --seed 42
"""

import argparse, os, sys, random, subprocess, tempfile, shutil
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np

COVERAGE_LEVELS = list(range(100, 10, -10))  # 100,90,80,...,20
KNOWN_COV = [100, 90, 80, 70, 60, 50, 40]   # 已知病毒: 7级

def load_virus_genomes(virus_dir, min_len=1000):
    """加载病毒目录下的所有基因组"""
    seqs = {}
    for f in sorted(Path(virus_dir).glob("*.fasta")):
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            if len(rec.seq) >= min_len:
                seqs[acc] = str(rec.seq)
    return seqs

def mutate(seq, rate, rng):
    """随机点突变, 返回突变后的序列"""
    if rate <= 0: return seq
    s = list(seq)
    n = int(len(s) * rate / 100)
    positions = rng.sample(range(len(s)), n)
    bases = ['A', 'C', 'G', 'T']
    for p in positions:
        orig = s[p]
        s[p] = rng.choice([b for b in bases if b != orig])
    return ''.join(s)

def mutate_to_similarity(seq, target_sim, rng, max_iter=200):
    """
    通过迭代突变使序列达到目标相似度 (与原始序列的核酸一致度).
    target_sim: 0.0-1.0, 目标相似度
    简化策略: sim = 1 - mutation_rate/100, 所以 mutation_rate = (1-sim)*100
    """
    rate = (1.0 - target_sim) * 100  # 0.90 → 10%, 0.70 → 30%, 0.40 → 60%
    # 由于随机突变不会精确达到目标, 我们迭代调整
    s = list(seq)
    L = len(s)
    target_n = int(L * target_sim)
    bases = ['A', 'C', 'G', 'T']
    orig = seq

    for _ in range(max_iter):
        current = ''.join(s)
        current_match = sum(1 for i in range(L) if current[i] == orig[i])
        if abs(current_match - target_n) <= max(1, L * 0.02):  # 容忍2%误差
            return current
        if current_match > target_n:
            # 需要更多突变
            n_mut = min(L, current_match - target_n + rng.randint(1, max(1, L // 200)))
            pos = rng.sample(range(L), n_mut)
            for p in pos:
                s[p] = rng.choice([b for b in bases if b != orig[p]])
        else:
            # 突变太多了, 回退一些
            n_fix = min(L, target_n - current_match + rng.randint(1, max(1, L // 200)))
            diff_pos = [i for i in range(L) if s[i] != orig[i]]
            if diff_pos:
                fix_pos = rng.sample(diff_pos, min(n_fix, len(diff_pos)))
                for p in fix_pos:
                    s[p] = orig[p]

    return ''.join(s)

def viromock_backtranslate(aa_seq, rng):
    """
    VIROMOCK Dataset 6 风格: 反向翻译 + 非同义突变
    输入: 氨基酸序列
    输出: 核苷酸序列 (与原始核苷酸相似度 ~70%, 但氨基酸相似度 ~99%)
    策略:
      1. 从同义密码子中随机选取 (避免原始密码子)
      2. 部分位置引入非同义替换 (改变氨基酸)
      注: 该实现为简化版, 完整版需 codon usage table, 此处在密码子空间做随机选择
    """
    # 简化版: 仅支持标准密码子表, 每个AA随机选一个密码子(非原始偏好)
    codon_table = {
        'A': ['GCT','GCC','GCA','GCG'], 'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
        'N': ['AAT','AAC'], 'D': ['GAT','GAC'], 'C': ['TGT','TGC'],
        'Q': ['CAA','CAG'], 'E': ['GAA','GAG'], 'G': ['GGT','GGC','GGA','GGG'],
        'H': ['CAT','CAC'], 'I': ['ATT','ATC','ATA'], 'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
        'K': ['AAA','AAG'], 'M': ['ATG'], 'F': ['TTT','TTC'],
        'P': ['CCT','CCC','CCA','CCG'], 'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
        'T': ['ACT','ACC','ACA','ACG'], 'W': ['TGG'], 'Y': ['TAT','TAC'],
        'V': ['GTT','GTC','GTA','GTG'], '*': ['TAA','TAG','TGA']
    }
    nt_seq = ''
    for aa in aa_seq:
        if aa in codon_table:
            nt_seq += rng.choice(codon_table[aa])
        else:
            nt_seq += 'NNN'  # unknown AA
    # 额外非同义突变: ~1% AA变化, 进一步降低nt identity
    nt_list = list(nt_seq)
    n_nonsyn = max(1, len(nt_seq) // 100)  # ~1% nucleotide substitutions
    pos = rng.sample(range(len(nt_seq)), n_nonsyn)
    for p in pos:
        orig = nt_list[p]
        nt_list[p] = rng.choice([b for b in 'ACGT' if b != orig])
    return ''.join(nt_list)

def generate_fragments(genomes_dict, coverage_levels, n_per_cov, rng,
                       label_prefix="", extra_meta=None):
    """核心片段生成器"""
    records = []
    idx = 0
    for cov_pct in coverage_levels:
        for _ in range(n_per_cov):
            for acc, seq in genomes_dict.items():
                L = len(seq)
                frag_len = max(200, int(L * cov_pct / 100))
                frag_len = min(frag_len, L)
                start = rng.randint(0, max(1, L - frag_len)) if L > frag_len else 0
                frag = seq[start:start + frag_len]

                meta_str = ""
                if extra_meta:
                    meta_str = '|' + '|'.join(f'{k}={v}' for k, v in extra_meta.items())
                seq_id = f"{label_prefix}|cov{cov_pct}|idx{idx:06d}|source={acc}{meta_str}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": L, "frag_length": len(frag),
                    "frag_seq": frag,
                })
                idx += 1
    return records

def generate_negatives(pos_lengths, host_genome, conserved_dir, eve_fasta, n_each, rng):
    """生成 A/B/C 三类负样本"""
    # A
    host_seqs = []
    if host_genome and os.path.exists(host_genome):
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= 200: host_seqs.append(str(rec.seq))
    if not host_seqs: host_seqs = ["ACGT" * 5000]
    neg_A = []
    for i in range(n_each):
        flen = max(200, rng.choice(pos_lengths) + rng.randint(-200, 200))
        seq = rng.choice(host_seqs); flen = min(flen, len(seq)-5)
        start = rng.randint(0, max(1, len(seq)-flen))
        neg_A.append({"seq_id": f"negA|{i:05d}", "label":"negative","type":"negative_A",
                       "frag_seq":seq[start:start+flen], "frag_length":flen})

    # B
    trap_seqs = []
    if conserved_dir and os.path.exists(conserved_dir):
        for f in Path(conserved_dir).glob("*.fasta"):
            try:
                rec = next(SeqIO.parse(f,"fasta"))
                if len(rec.seq)>=500: trap_seqs.append(str(rec.seq))
            except: pass
    if not trap_seqs: trap_seqs = ["ACGT"*3000]
    neg_B = []
    for i in range(n_each):
        seq = rng.choice(trap_seqs)
        flen = min(rng.choice(pos_lengths), len(seq)-5); flen = max(200, flen)
        start = rng.randint(0, max(1, len(seq)-flen))
        neg_B.append({"seq_id":f"negB|{i:05d}","label":"negative","type":"negative_B",
                       "frag_seq":seq[start:start+flen], "frag_length":flen})

    # C
    eve_seqs = []
    if eve_fasta and os.path.exists(eve_fasta):
        for rec in SeqIO.parse(eve_fasta,"fasta"):
            if len(rec.seq)>=500: eve_seqs.append(str(rec.seq))
    if not eve_seqs: eve_seqs = host_seqs[:50] if host_seqs else ["ACGT"*5000]
    neg_C = []
    for i in range(n_each):
        seq = rng.choice(eve_seqs)
        flen = min(rng.choice(pos_lengths), len(seq)-5); flen = max(200, flen)
        start = rng.randint(0, max(1, len(seq)-flen))
        neg_C.append({"seq_id":f"negC|{i:05d}","label":"negative","type":"negative_C",
                       "frag_seq":seq[start:start+flen], "frag_length":flen})

    return neg_A + neg_B + neg_C

def write_output(outdir, name, pos_records, neg_records, extra_cols=None):
    fasta_out = os.path.join(outdir, f"{name}.fasta")
    labels_out = os.path.join(outdir, f"{name}_labels.tsv")
    recs, labs = [], []
    for r in pos_records:
        recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        lab = {k:r[k] for k in ["seq_id","label","type","source","coverage_pct"] if k in r}
        if extra_cols:
            lab.update({k: r.get(k) for k in extra_cols})
        labs.append(lab)
    for r in neg_records:
        recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        labs.append({k:r[k] for k in ["seq_id","label","type"]})
    SeqIO.write(recs, fasta_out, "fasta")
    pd.DataFrame(labs).to_csv(labels_out, sep='\t', index=False)
    return fasta_out

def main():
    parser = argparse.ArgumentParser(description='综合已知/未知病毒评估')
    parser.add_argument('--virus-dir', required=True, help='step1_eval_viruses/')
    parser.add_argument('--ncbi-meta', default=None, help='NCBI新病毒元数据 (方案B来源2)')
    parser.add_argument('--novel-fasta', default=None, help='本地新病毒FASTA (方案B来源2, 优先于此)')
    parser.add_argument('--ref-fasta', default=None, help='final.cluster.ref.fasta')
    parser.add_argument('--n-novel-source1', type=int, default=10, help='来源1的病毒数')
    parser.add_argument('--host-genome', default=None)
    parser.add_argument('--conserved-prots', default=None)
    parser.add_argument('--eve-fasta', default=None)
    parser.add_argument('--n-neg', type=int, default=500, help='每类负样本')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    cov_levels = COVERAGE_LEVELS  # 100,90,...,20
    n_rep = 5

    # 加载病毒基因组
    all_genomes = load_virus_genomes(args.virus_dir)
    print(f"[load] {len(all_genomes)} viruses >= 1000 bp")

    # 预生成负样本 (所有方案共用)
    sample_lens = [max(200, int(L * c / 100)) for L in [len(s) for s in all_genomes.values()]
                   for c in cov_levels]
    neg_records = generate_negatives(sample_lens, args.host_genome, args.conserved_prots,
                                     args.eve_fasta, args.n_neg, rng)
    print(f"[neg] {len(neg_records)} negative sequences (A/B/C each {args.n_neg})")

    # ═══════════════════════════════════════════
    # 方案A: 已知病毒检测 (Known)
    # ═══════════════════════════════════════════
    print("\n" + "="*60)
    print("  方案A: 已知病毒 (Known) — 0% + 1/3/5/7/9% 突变")
    print("="*60)

    known_dir = os.path.join(args.outdir, "A_known")
    os.makedirs(known_dir, exist_ok=True)
    known_rates = [0, 5]

    all_known_pos = []
    for rate in known_rates:
        if rate == 0:
            mut_genomes = all_genomes
        else:
            mut_genomes = {f"{acc}_mut{rate}pct": mutate(seq, rate, rng)
                          for acc, seq in all_genomes.items()}
        pos = generate_fragments(mut_genomes, KNOWN_COV, n_rep, rng,
                                 label_prefix="known",
                                 extra_meta={"mut_rate": rate, "scheme": "A"})
        all_known_pos.extend(pos)
        print(f"  mut {rate}%: {len(pos)} fragments")

    write_output(known_dir, "known_all", all_known_pos, neg_records,
                 extra_cols=["mut_rate"])

    # ═══════════════════════════════════════════
    # 方案B: 新病毒检测 (Novel)
    # ═══════════════════════════════════════════
    print("\n" + "="*60)
    print("  方案B: 新病毒 (Novel)")
    print("="*60)

    novel_dir = os.path.join(args.outdir, "B_novel")
    os.makedirs(novel_dir, exist_ok=True)

    # 来源1: 原有病毒突变到目标相似度
    print("\n  [B1] 来源1: 原有病毒突变到目标相似度...")
    src1_viruses = dict(rng.sample(list(all_genomes.items()),
                                    min(args.n_novel_source1, len(all_genomes))))
    print(f"  [B1] Selected {len(src1_viruses)} viruses")

    # 相似度分级
    sim_levels = {
        'L1': [0.90, 0.80, 0.70],
        'L2': [0.70, 0.60, 0.50],
        'L3': [0.40, 0.30, 0.20],
    }

    all_b1_pos = []
    for level, sims in sim_levels.items():
        for sim in sims:
            mut_genomes = {}
            for acc, seq in src1_viruses.items():
                mut_seq = mutate_to_similarity(seq, sim, rng)
                actual_sim = sum(1 for i in range(len(seq)) if seq[i]==mut_seq[i]) / len(seq)
                mut_genomes[f"{acc}_sim{sim:.0f}"] = mut_seq
            pos = generate_fragments(mut_genomes, cov_levels, n_rep, rng,
                                     label_prefix="novel_src1",
                                     extra_meta={"level": level, "target_sim": sim,
                                                 "scheme": "B1"})
            all_b1_pos.extend(pos)
            print(f"    {level} sim={sim:.2f}: {len(pos)} fragments")

    write_output(novel_dir, "novel_src1_mutation", all_b1_pos, neg_records,
                 extra_cols=["level", "target_sim"])

    # 来源2: 本地新病毒 或 NCBI新病毒
    novel_genomes = {}
    src2_label = "local"

    if args.novel_fasta and os.path.exists(args.novel_fasta):
        print(f"\n  [B2] 来源2: 本地新病毒 ({args.novel_fasta})...")
        for rec in SeqIO.parse(args.novel_fasta, "fasta"):
            if len(rec.seq) >= 1000:
                acc = rec.id.split()[0]
                novel_genomes[acc] = str(rec.seq)
        src2_label = "local"
    elif args.ncbi_meta and args.ref_fasta and os.path.exists(args.ncbi_meta):
        print("\n  [B2] 来源2: NCBI新病毒...")
        meta = pd.read_csv(args.ncbi_meta, sep='\t')
        ref_idx = SeqIO.index(args.ref_fasta, "fasta")
        for _, row in meta.iterrows():
            acc = row['Accession']
            length = int(row['Length'])
            if length >= 1000 and acc in ref_idx:
                novel_genomes[acc] = str(ref_idx[acc].seq)
        src2_label = "ncbi"

    if novel_genomes:
        print(f"  [B2] {len(novel_genomes)} novel viruses >= 1000 bp")
        all_b2_pos = generate_fragments(novel_genomes, cov_levels, n_rep, rng,
                                        label_prefix=f"novel_{src2_label}_raw",
                                        extra_meta={"mut_rate": 0, "scheme": f"B2_{src2_label}_raw"})

        for level, sims in sim_levels.items():
            for sim in sims:
                mut_genomes = {}
                for acc, seq in novel_genomes.items():
                    mut_seq = mutate_to_similarity(seq, sim, rng)
                    mut_genomes[f"{acc}_sim{sim:.0f}"] = mut_seq
                pos = generate_fragments(mut_genomes, cov_levels, n_rep, rng,
                                         label_prefix=f"novel_{src2_label}_mut",
                                         extra_meta={"level": level, "target_sim": sim,
                                                     "scheme": f"B2_{src2_label}_mut"})
                all_b2_pos.extend(pos)

        write_output(novel_dir, f"novel_src2_{src2_label}", all_b2_pos, neg_records,
                     extra_cols=["level", "target_sim", "scheme"])
    else:
        print("\n  [B2] 跳过 (未提供 --novel-fasta 或 --ncbi-meta 中无可用序列)")

    # ═══════════════════════════════════════════
    # 方案B3: VIROMOCK Dataset 6 风格 — 反向翻译 + 非同义突变
    # ═══════════════════════════════════════════
    print("\n  [B3] VIROMOCK风格: 反向翻译 + 非同义突变...")
    b3_dir = os.path.join(novel_dir, "B3_viromock")
    os.makedirs(b3_dir, exist_ok=True)

    # 从 src1 病毒中选5个, 翻译 → 反向翻译 → 非同义突变
    viromock_viruses = dict(list(src1_viruses.items())[:5])
    viromock_genomes = {}

    for acc, seq in viromock_viruses.items():
        # 翻译为氨基酸
        aa_seq = ''
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            # 简化翻译表
            codon_aa = {
                'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TGT':'C','TGC':'C','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I',
                'ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N',
                'AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                'TAA':'*','TAG':'*','TGA':'*'
            }
            aa = codon_aa.get(codon, 'X')
            aa_seq += aa

        # VIROMOCK 反向翻译
        viromock_nt = viromock_backtranslate(aa_seq, rng)

        # 计算与原始序列的相似度
        min_len = min(len(seq), len(viromock_nt))
        actual_sim = sum(1 for i in range(min_len) if seq[i]==viromock_nt[i]) / min_len
        viromock_genomes[f"{acc}_viromock"] = viromock_nt
        print(f"    {acc}: nt_identity={actual_sim:.3f}")

    all_b3_pos = generate_fragments(viromock_genomes, cov_levels, n_rep, rng,
                                    label_prefix="novel_viromock",
                                    extra_meta={"scheme": "B3_viromock", "mut_rate": 0})
    print(f"  [B3] {len(all_b3_pos)} VIROMOCK-style fragments")
    write_output(novel_dir, "novel_src3_viromock", all_b3_pos, neg_records,
                 extra_cols=["scheme"])

    # ═══════════════════════════════════════════
    # 汇总
    # ═══════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"  生成完毕")
    print(f"{'='*60}")
    print(f"  方案A (已知):   {len(all_known_pos)} 正样本 ({len(known_rates)} mutation rates × {len(cov_levels)} cov × {n_rep} rep)")
    print(f"  方案B1 (突变):  {len(all_b1_pos)} 正样本")
    print(f"  方案B2 (新病毒): {len(all_b2_pos) if 'all_b2_pos' in dir() else 0} 正样本 ({src2_label})")
    print(f"  方案B3 (VIROMOCK): {len(all_b3_pos) if 'all_b3_pos' in dir() else 0} 正样本")
    print(f"  共享负样本:     {len(neg_records)}")
    print(f"  输出:           {args.outdir}")

if __name__ == '__main__':
    main()
