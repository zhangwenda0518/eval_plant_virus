#!/usr/bin/env python3
"""
统一模拟数据集构建器

方案:
  A  已知病毒:  已知60+新45 → 突变 0%/5%/10% → 9覆盖率 → 5重复
  B  新病毒(突变): 已知60+新45 → 相似度 80/70/60/50/40/30/20 → 9覆盖率 → 5重复
  C  新病毒(VIROMOCK): 已知60+新45 → 反向翻译 → 相似度 80/70/60/50/40/30/20 → 9覆盖率 → 5重复
  D  Conserved_Traps:  1000
  E  EVE_Transposon:   1000
  H  Host_Random:      1000

标签命名:
  A: A{mut}C{cov}R{rep}  例: A0C100R1, A5C90R2, A10C40R5
  B: B{sim}C{cov}R{rep}  例: B80C100R1, B50C60R3
  C: C{sim}C{cov}R{rep}  例: C80C100R1, C50C60R3
  D/E/H: D/E/H + 序号

用法:
  python prep_simulation_datasets.py \
      --known-viruses step1_eval_viruses/ \
      --new-viruses new_viruses.fasta \
      --host-genome ningxia.genome.fasta \
      --conserved-prots step3_conserved_traps/ \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --outdir simulation_datasets/ --seed 42
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
KNOWN_MUT_RATES = [0, 5, 10]
NOVEL_SIM_LEVELS = [80, 70, 60, 50, 40, 30, 20]
MIN_LEN = 1000

# 密码子表
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
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def load_genomes(path):
    """加载FASTA目录或文件, 返回 {acc: seq}"""
    genomes = {}
    loader = [Path(path).glob("*.fasta")] if os.path.isdir(path) else [[Path(path)]]
    for files in loader:
        for f in sorted(files):
            try:
                for rec in SeqIO.parse(f, "fasta"):
                    if len(rec.seq) >= MIN_LEN:
                        genomes[rec.id.split()[0]] = str(rec.seq)
            except: pass
    return genomes


def mutate_random(seq, rate_pct, rng):
    """随机点突变 rate_pct%"""
    if rate_pct <= 0: return seq
    s = list(seq)
    n = int(len(s) * rate_pct / 100)
    for p in rng.sample(range(len(s)), min(n, len(s))):
        orig = s[p]
        s[p] = rng.choice([b for b in 'ACGT' if b != orig])
    return ''.join(s)


def mutate_to_similarity(seq, target_pct, rng, max_iter=200):
    """迭代突变达到目标相似度 target_pct%"""
    target_match = int(len(seq) * target_pct / 100)
    s = list(seq)
    for _ in range(max_iter):
        current_match = sum(1 for i in range(len(seq)) if s[i] == seq[i])
        if abs(current_match - target_match) <= max(1, len(seq) * 0.02):
            return ''.join(s)
        if current_match > target_match:
            n = min(len(seq), current_match - target_match + 1)
            for p in rng.sample(range(len(seq)), n):
                s[p] = rng.choice([b for b in 'ACGT' if b != seq[p]])
        else:
            diff = [i for i in range(len(seq)) if s[i] != seq[i]]
            if diff:
                n = min(len(diff), target_match - current_match + 1)
                for p in rng.sample(diff, n):
                    s[p] = seq[p]
    return ''.join(s)


def translate_to_aa(seq):
    """翻译为氨基酸 (frame 0, 忽略提前终止)"""
    aa = ''
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa += CODON_AA.get(codon, 'X')
    return aa


def backtranslate_to_similarity(orig_seq, target_pct, rng):
    """
    VIROMOCK 反向翻译到目标相似度
    编码区: 密码子重选 + 安全点突变
    非编码区: 简单点突变
    简化版: 将整个序列视为一个大CDS
    """
    S_base = 68  # 纯随机同义反向翻译基线 ~68%

    if target_pct > S_base:
        p_keep = (target_pct - S_base) / (100 - S_base)
        extra_mut = 0.0
    else:
        p_keep = 0.0
        extra_mut = (S_base - target_pct) / S_base

    codons = [orig_seq[i:i+3] for i in range(0, len(orig_seq) - 2, 3)]
    new_codons = []
    for codon in codons:
        if codon not in CODON_AA:
            new_codons.append(codon)
            continue
        if rng.random() < p_keep:
            new_codons.append(codon)
        else:
            aa = CODON_AA[codon]
            choices = list(CODON_USAGE[aa].keys())
            weights = list(CODON_USAGE[aa].values())
            new_codons.append(rng.choices(choices, weights=weights, k=1)[0])

    result = ''.join(new_codons)
    if extra_mut > 0:
        result = apply_safe_mutations(result, extra_mut, rng)
    return result


def apply_safe_mutations(seq, rate, rng):
    """编码区安全突变: 避免创建提前终止密码子"""
    s = list(seq)
    n = int(len(s) * rate)
    for p in rng.sample(range(len(s)), min(n, len(s))):
        codon_start = (p // 3) * 3
        if codon_start + 3 > len(s): continue
        old_codon = ''.join(s[codon_start:codon_start+3])
        p_in = p % 3
        valid = []
        for b in 'ACGT':
            if b == s[p]: continue
            new_codon = list(old_codon); new_codon[p_in] = b
            if old_codon in STOP_CODONS:
                if ''.join(new_codon) in STOP_CODONS: valid.append(b)
            else:
                if ''.join(new_codon) not in STOP_CODONS: valid.append(b)
        if valid: s[p] = rng.choice(valid)
    return ''.join(s)


def generate_fragments(genomes, cov_levels, n_rep, rng, prefix, extra_meta_func=None):
    """生成片段: prefix=标签前缀, extra_meta_func=lambda acc: 额外元信息dict"""
    records = []
    idx = 0
    for cov in cov_levels:
        for rep in range(1, n_rep + 1):
            for acc, seq in genomes.items():
                L = len(seq)
                flen = max(200, int(L * cov / 100))
                flen = min(flen, L)
                start = rng.randint(0, L - flen) if L > flen else 0
                frag = seq[start:start + flen]

                extra = extra_meta_func(acc) if extra_meta_func else {}
                meta = '|'.join(f'{k}={v}' for k, v in extra.items())
                meta_str = f'|{meta}' if meta else ''
                seq_id = f"{prefix}C{cov}R{rep}|idx{idx:06d}|source={acc}{meta_str}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov, "replicate": rep,
                    "frag_seq": frag, "frag_length": flen,
                })
                idx += 1
    return records


def generate_negatives(pos_lengths_pool, host_genome, conserved_dir, eve_fasta, n_each, rng):
    """生成 D/E/H 三类负样本"""
    results = {}

    # H: Host_Random
    host = [] if not host_genome or not os.path.exists(host_genome) else [
        str(r.seq) for r in SeqIO.parse(host_genome, "fasta") if len(r.seq) >= 200]
    if not host: host = ["ACGT" * 5000]
    negs_h = []
    for i in range(n_each):
        tlen = rng.choice(pos_lengths_pool)
        flen = max(200, tlen + rng.randint(-200, 200))
        s = rng.choice(host); flen = min(flen, len(s)-5)
        st = rng.randint(0, len(s)-flen)
        negs_h.append({"seq_id": f"H{chr(65+i//1000)}{i%1000:03d}", "label":"negative",
                       "type":"Host_Random", "frag_seq":s[st:st+flen], "frag_length":flen})
    results['H'] = negs_h

    # D: Conserved_Traps
    traps = []
    if conserved_dir and os.path.exists(conserved_dir):
        for f in Path(conserved_dir).glob("*.fasta"):
            try:
                r = next(SeqIO.parse(f,"fasta"))
                if len(r.seq)>=500: traps.append(str(r.seq))
            except: pass
    if not traps: traps = ["ACGT"*3000]
    negs_d = []
    for i in range(n_each):
        s = rng.choice(traps)
        flen = min(rng.choice(pos_lengths_pool), len(s)-5)
        flen = max(200, flen)
        st = rng.randint(0, len(s)-flen)
        negs_d.append({"seq_id": f"D{chr(65+i//1000)}{i%1000:03d}", "label":"negative",
                       "type":"Conserved_Traps", "frag_seq":s[st:st+flen], "frag_length":flen})
    results['D'] = negs_d

    # E: EVE_Transposon
    eves = []
    if eve_fasta and os.path.exists(eve_fasta):
        for r in SeqIO.parse(eve_fasta,"fasta"):
            if len(r.seq)>=500: eves.append(str(r.seq))
    if not eves: eves = host[:50] if host else ["ACGT"*5000]
    negs_e = []
    for i in range(n_each):
        s = rng.choice(eves)
        flen = min(rng.choice(pos_lengths_pool), len(s)-5)
        flen = max(200, flen)
        st = rng.randint(0, len(s)-flen)
        negs_e.append({"seq_id": f"E{chr(65+i//1000)}{i%1000:03d}", "label":"negative",
                       "type":"EVE_Transposon", "frag_seq":s[st:st+flen], "frag_length":flen})
    results['E'] = negs_e

    return results


def merge_and_write(outdir, name, all_pos, neg_sets, extra_cols=None):
    """合并正负样本写入FASTA+labels"""
    recs, labs = [], []
    for r in all_pos:
        recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        lab = {k: r[k] for k in ["seq_id","label","type","source","coverage_pct"] if k in r}
        if extra_cols:
            lab.update({k: r.get(k) for k in extra_cols})
        labs.append(lab)
    for neg_type, negs in neg_sets.items():
        for r in negs:
            recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
            labs.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})

    SeqIO.write(recs, os.path.join(outdir, f"{name}.fasta"), "fasta")
    pd.DataFrame(labs).to_csv(os.path.join(outdir, f"{name}_labels.tsv"), sep='\t', index=False)
    return len(recs)


def main():
    parser = argparse.ArgumentParser(description='统一模拟评估数据集构建器')
    parser.add_argument('--known-viruses', required=True, help='已知60病毒的FASTA目录')
    parser.add_argument('--new-viruses', required=True, help='新45病毒的FASTA文件')
    parser.add_argument('--host-genome', default=None)
    parser.add_argument('--conserved-prots', default=None)
    parser.add_argument('--eve-fasta', default=None)
    parser.add_argument('--n-neg', type=int, default=1000, help='每类负样本数')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    t0 = time.time()
    os.makedirs(args.outdir, exist_ok=True)

    # ── 加载基因组 ──
    known = load_genomes(args.known_viruses)
    novel = load_genomes(args.new_viruses)
    all_virus = {**known, **novel}
    print(f"[load] Known: {len(known)}  Novel: {len(novel)}  Total: {len(all_virus)}")

    # ── 预生成负样本(共用) ──
    pos_lens_pool = [max(200, int(L*c/100)) for L in [len(s) for s in all_virus.values()]
                     for c in COV_LEVELS]
    neg_sets = generate_negatives(pos_lens_pool, args.host_genome,
                                  args.conserved_prots, args.eve_fasta, args.n_neg, rng)
    for k, v in neg_sets.items():
        print(f"[neg] {k}: {len(v)}")

    # ═══════════════════════════════════════════
    # 方案A: 已知病毒评估
    # ═══════════════════════════════════════════
    print(f"\n{'='*60}\n  方案A: 已知病毒 (突变 0%/5%/10%)\n{'='*60}")
    pos_a = []
    for rate in KNOWN_MUT_RATES:
        if rate == 0:
            mut_genomes = all_virus
        else:
            mut_genomes = {f"{acc}_mut{rate}pct": mutate_random(seq, rate, rng)
                           for acc, seq in all_virus.items()}
        frags = generate_fragments(mut_genomes, COV_LEVELS, N_REP, rng,
                                   prefix=f"A{rate}",
                                   extra_meta_func=lambda a, r=rate: {"mut_rate": r, "scheme": "A"})
        pos_a.extend(frags)
        print(f"  A mut={rate}%: {len(frags)} fragments")

    n_a = merge_and_write(args.outdir, "A_known", pos_a, neg_sets, ["mut_rate", "scheme"])

    # ═══════════════════════════════════════════
    # 方案B: 新病毒(突变相似度)
    # ═══════════════════════════════════════════
    print(f"\n{'='*60}\n  方案B: 新病毒突变 (相似度 80-20%)\n{'='*60}")
    pos_b = []
    for sim in NOVEL_SIM_LEVELS:
        # 对每个病毒, 预计算突变(相同sim共享)
        sim_genomes = {}
        for acc, seq in all_virus.items():
            sim_genomes[f"{acc}_sim{sim}"] = mutate_to_similarity(seq, sim, rng)
        frags = generate_fragments(sim_genomes, COV_LEVELS, N_REP, rng,
                                   prefix=f"B{sim}",
                                   extra_meta_func=lambda a, s=sim: {"target_sim": s, "scheme": "B"})
        pos_b.extend(frags)
        # 验证一个实际相似度
        sample_acc = list(sim_genomes.keys())[0]
        sample_seq = sim_genomes[sample_acc]
        orig_acc = sample_acc.rsplit('_sim', 1)[0]
        actual = sum(1 for i in range(min(len(sample_seq), len(all_virus[orig_acc])))
                     if sample_seq[i] == all_virus[orig_acc][i]) / min(len(sample_seq), len(all_virus[orig_acc]))
        print(f"  B sim={sim}%: {len(frags)} fragments (actual_sim~{actual*100:.0f}%)")

    n_b = merge_and_write(args.outdir, "B_novel_mutation", pos_b, neg_sets, ["target_sim", "scheme"])

    # ═══════════════════════════════════════════
    # 方案C: VIROMOCK 反向翻译
    # ═══════════════════════════════════════════
    print(f"\n{'='*60}\n  方案C: VIROMOCK 反向翻译 (相似度 80-20%)\n{'='*60}")
    pos_c = []
    for sim in NOVEL_SIM_LEVELS:
        vm_genomes = {}
        for acc, seq in all_virus.items():
            vm_seq = backtranslate_to_similarity(seq, sim, rng)
            vm_genomes[f"{acc}_vm{sim}"] = vm_seq
        frags = generate_fragments(vm_genomes, COV_LEVELS, N_REP, rng,
                                   prefix=f"C{sim}",
                                   extra_meta_func=lambda a, s=sim: {"target_sim": s, "scheme": "C"})
        pos_c.extend(frags)
        sample_acc = list(vm_genomes.keys())[0]
        sample_seq = vm_genomes[sample_acc]
        orig_acc = sample_acc.rsplit('_vm', 1)[0]
        actual = sum(1 for i in range(min(len(sample_seq), len(all_virus[orig_acc])))
                     if sample_seq[i] == all_virus[orig_acc][i]) / min(len(sample_seq), len(all_virus[orig_acc]))
        print(f"  C sim={sim}%: {len(frags)} fragments (actual_sim~{actual*100:.0f}%)")

    n_c = merge_and_write(args.outdir, "C_viromock", pos_c, neg_sets, ["target_sim", "scheme"])

    # ═══════════════════════════════════════════
    # 汇总
    # ═══════════════════════════════════════════
    total_pos = len(pos_a) + len(pos_b) + len(pos_c)
    total_neg = sum(len(v) for v in neg_sets.values())
    elapsed = (time.time() - t0) / 60
    print(f"\n{'='*60}")
    print(f"  生成完毕 ({elapsed:.1f} min)")
    print(f"{'='*60}")
    print(f"  A 已知病毒:         {len(pos_a):>8,} 正样本 ({len(all_virus)} viruses × {len(KNOWN_MUT_RATES)} rates × {len(COV_LEVELS)} cov × {N_REP} rep)")
    print(f"  B 新病毒(突变):     {len(pos_b):>8,} 正样本 ({len(all_virus)} viruses × {len(NOVEL_SIM_LEVELS)} sim × {len(COV_LEVELS)} cov × {N_REP} rep)")
    print(f"  C 新病毒(VIROMOCK): {len(pos_c):>8,} 正样本 ({len(all_virus)} viruses × {len(NOVEL_SIM_LEVELS)} sim × {len(COV_LEVELS)} cov × {N_REP} rep)")
    print(f"  正样本总计:         {total_pos:>8,}")
    print(f"  负样本:             {total_neg:>8,}")
    print(f"  全部序列:           {total_pos+total_neg:>8,}")
    print(f"  输出:               {args.outdir}")

if __name__ == '__main__':
    main()
