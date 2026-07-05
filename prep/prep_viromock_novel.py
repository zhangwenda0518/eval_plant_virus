#!/usr/bin/env python3
"""
VIROMOCK Dataset 6 风格的新病毒序列生成器 (独立脚本)
参考: Tamisier et al. (2021) Peer Community Journal, 1: e53

原理:
  1. 将病毒基因组翻译为氨基酸序列
  2. 反向翻译: 从同义密码子中随机选取 (利用密码子简并性)
     → 氨基酸100%一致, 核苷酸~70%一致
  3. 额外引入~1%非同义突变
     → 氨基酸~99%一致, 核苷酸~70%一致
  4. 按覆盖率梯度截取片段, 输出评估用FASTA

用法:
  # 第一步: 准备输入病毒FASTA (任意来源, 越多越好)
  # 第二步: 运行
  python prep_viromock_novel.py \
      --input step1_eval_viruses/ \
      --n-viruses 10 \
      --coverage-levels 100 80 60 40 \
      --n-per-coverage 5 \
      --outdir step3_viromock/ --seed 42

输出:
  {outdir}/
    ├── viromock_genomes.fasta    # 反向翻译后的全长基因组
    ├── viromock_fragments.fasta   # 覆盖率梯度片段
    └── viromock_labels.tsv        # 标签文件
"""

import argparse, os, sys, random
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np


# 标准遗传密码子表
CODON_TABLE = {
    'A': ['GCT','GCC','GCA','GCG'], 'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
    'N': ['AAT','AAC'], 'D': ['GAT','GAC'], 'C': ['TGT','TGC'],
    'Q': ['CAA','CAG'], 'E': ['GAA','GAG'], 'G': ['GGT','GGC','GGA','GGG'],
    'H': ['CAT','CAC'], 'I': ['ATT','ATC','ATA'], 'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
    'K': ['AAA','AAG'], 'M': ['ATG'], 'F': ['TTT','TTC'],
    'P': ['CCT','CCC','CCA','CCG'], 'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
    'T': ['ACT','ACC','ACA','ACG'], 'W': ['TGG'], 'Y': ['TAT','TAC'],
    'V': ['GTT','GTC','GTA','GTG'], '*': ['TAA','TAG','TGA']
}

CODON_AA = {
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


def translate(seq):
    """将核苷酸序列翻译为氨基酸序列 (正向3框, 取最长ORF)"""
    best_aa = ''
    for frame in range(3):
        aa = ''
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            aa += CODON_AA.get(codon, 'X')
        # 取第一个 > 50% 非X的
        if sum(1 for c in aa if c != 'X') > len(aa) * 0.5:
            return aa
    # fallback: frame 0
    aa = ''
    for i in range(0, len(seq) - 2, 3):
        aa += CODON_AA.get(seq[i:i+3], 'X')
    return aa


def backtranslate(aa_seq, rng):
    """
    VIROMOCK反向翻译: 每个氨基酸随机选取同义密码子
    结果: 氨基酸序列100%一致, 核苷酸序列~70%一致
    """
    nt_seq = ''
    for aa in aa_seq:
        if aa in CODON_TABLE:
            nt_seq += rng.choice(CODON_TABLE[aa])
        else:
            nt_seq += 'NNN'
    return nt_seq


def add_nonsynonymous(nt_seq, n_substitutions, rng):
    """随机引入非同义替换 (~1% 氨基酸变化)"""
    nt_list = list(nt_seq)
    positions = rng.sample(range(len(nt_seq)), min(n_substitutions, len(nt_seq)))
    for p in positions:
        orig = nt_list[p]
        nt_list[p] = rng.choice([b for b in 'ACGT' if b != orig])
    return ''.join(nt_list)


def generate_fragments(genomes, coverage_levels, n_per_cov, rng):
    """按覆盖率梯度截取片段"""
    records = []
    idx = 0
    for cov_pct in coverage_levels:
        for _ in range(n_per_cov):
            for acc, seq in genomes.items():
                L = len(seq)
                frag_len = max(200, int(L * cov_pct / 100))
                frag_len = min(frag_len, L)
                start = rng.randint(0, max(1, L - frag_len)) if L > frag_len else 0
                frag = seq[start:start + frag_len]
                seq_id = f"viromock|cov{cov_pct}|idx{idx:05d}|source={acc}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": L, "frag_length": len(frag),
                    "frag_seq": frag,
                })
                idx += 1
    return records


def main():
    parser = argparse.ArgumentParser(description='VIROMOCK Dataset 6 风格新病毒生成器')
    parser.add_argument('--input', required=True,
                        help='输入病毒FASTA文件或目录 (≥1000 bp的病毒基因组)')
    parser.add_argument('--n-viruses', type=int, default=10,
                        help='选取的病毒数量 (从输入中随机选)')
    parser.add_argument('--coverage-levels', type=int, nargs='+',
                        default=[100, 80, 60, 40],
                        help='覆盖率等级 (默认: 100 80 60 40)')
    parser.add_argument('--n-per-coverage', type=int, default=5,
                        help='每个覆盖率等级的片段数 (默认: 5)')
    parser.add_argument('--nonsyn-rate', type=float, default=0.01,
                        help='非同义突变比例 (默认: 0.01 = 1%%)')
    parser.add_argument('--outdir', required=True, help='输出目录')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # 1. 加载输入病毒
    genomes = {}
    input_path = args.input
    if os.path.isdir(input_path):
        for f in sorted(Path(input_path).glob("*.fasta")):
            for rec in SeqIO.parse(f, "fasta"):
                if len(rec.seq) >= 1000:
                    genomes[rec.id.split()[0]] = str(rec.seq)
    elif os.path.isfile(input_path):
        for rec in SeqIO.parse(input_path, "fasta"):
            if len(rec.seq) >= 1000:
                genomes[rec.id.split()[0]] = str(rec.seq)

    print(f"[load] {len(genomes)} genomes >= 1000 bp")

    # 选取n个
    selected = dict(rng.sample(list(genomes.items()),
                                min(args.n_viruses, len(genomes))))
    print(f"[select] {len(selected)} viruses for VIROMOCK processing")

    # 2. 翻译 → 反向翻译 → 非同义突变
    viromock_genomes = {}
    viromock_records = []

    for acc, seq in selected.items():
        aa_seq = translate(seq)
        viromock_nt = backtranslate(aa_seq, rng)
        # 非同义突变
        n_sub = max(1, int(len(viromock_nt) * args.nonsyn_rate))
        viromock_nt = add_nonsynonymous(viromock_nt, n_sub, rng)

        # 计算与原始序列的相似度
        min_len = min(len(seq), len(viromock_nt))
        matches = sum(1 for i in range(min_len) if seq[i] == viromock_nt[i])
        nt_identity = matches / min_len if min_len > 0 else 0.0
        # 氨基酸相似度 (反向翻译前后)
        aa_orig = translate(seq)
        aa_new = translate(viromock_nt)
        min_aa = min(len(aa_orig), len(aa_new))
        aa_identity = sum(1 for i in range(min_aa) if aa_orig[i] == aa_new[i]) / min_aa if min_aa > 0 else 0.0

        viromock_genomes[acc] = viromock_nt
        viromock_records.append(SeqRecord(Seq(viromock_nt),
                                          id=f"{acc}_viromock",
                                          description=f"nt_id={nt_identity:.3f} aa_id={aa_identity:.3f}"))
        print(f"  {acc}: nt_id={nt_identity:.3f}  aa_id={aa_identity:.3f}  len={len(seq)}→{len(viromock_nt)}")

    # 4. 输出全长基因组
    genome_out = os.path.join(args.outdir, "viromock_genomes.fasta")
    SeqIO.write(viromock_records, genome_out, "fasta")
    print(f"\n[genomes] {genome_out}")

    # 5. 生成覆盖率梯度片段
    fragments = generate_fragments(viromock_genomes, args.coverage_levels,
                                   args.n_per_coverage, rng)
    print(f"[fragments] {len(fragments)} positive fragments "
          f"({len(selected)} viruses × {len(args.coverage_levels)} cov × {args.n_per_coverage} rep)")

    # 6. 输出
    recs = [SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description="")
            for r in fragments]
    frag_out = os.path.join(args.outdir, "viromock_fragments.fasta")
    SeqIO.write(recs, frag_out, "fasta")

    labels = [{"seq_id": r["seq_id"], "label": r["label"], "type": r["type"],
               "source": r["source"], "coverage_pct": r["coverage_pct"]}
              for r in fragments]
    labels_out = os.path.join(args.outdir, "viromock_labels.tsv")
    pd.DataFrame(labels).to_csv(labels_out, sep='\t', index=False)

    print(f"[output] {frag_out}")
    print(f"[output] {labels_out}")
    print(f"\n  下一步: 将 viromock_fragments.fasta 与负样本合并后提供给评估工具")


if __name__ == '__main__':
    main()
