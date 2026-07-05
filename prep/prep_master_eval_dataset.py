#!/usr/bin/env python3
"""
MASTER 版未知病毒鉴定评估数据集一键构建脚本 (本地序列注入升级版)

配置：
  1. 方案 A (已知): 已知病毒60条, 突变率 [0%, 5%, 10%], 覆盖度 [100, 90, 70, 50, 30, 10], 重复 5
  2. 方案 B (新病毒): 真实新病毒 (45条 + 本地追加序列), 0% 原始状态, 覆盖度 [100, 90, 70, 50, 30, 10], 重复 5
  3. 方案 C (VIROMOCK): 全部病毒, 70% 相似度, 严格 100% 全长, 5次独立密码子重组回译
  4. 负样本 (D/E/H): 各 1000 条 (默认值)

ID格式：
  - A/B方案: positive|scheme_[A|B]|mut[X]|cov[Y]|rep[Z]|source=[Accession]
  - C方案  : positive|scheme_C|sim70|cov100|rep[Z]|source=[Accession]
  - 负样本 : negative|[conserved_trap|EVE|host_random]|idx[Idx]|source=[Source]
"""

import argparse
import os
import sys
import random
import tempfile
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
try:
    import pyrodigal
    HAS_PYRODIGAL = True
except ImportError:
    HAS_PYRODIGAL = False
    print("[WARN] pyrodigal not installed, VIROMOCK will use simplified codon-reshuffling (no CDS segmentation)")

# 简并密码子偏好权重表
CODON_USAGE_BIAS = {
    'A': {'GCT': 0.30, 'GCC': 0.40, 'GCA': 0.15, 'GCG': 0.15},
    'R': {'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20, 'AGG': 0.21},
    'N': {'AAT': 0.46, 'AAC': 0.54}, 'D': {'GAT': 0.54, 'GAC': 0.46},
    'C': {'TGT': 0.45, 'TGC': 0.55}, 'Q': {'CAA': 0.25, 'CAG': 0.75},
    'E': {'GAA': 0.42, 'GAG': 0.58}, 'G': {'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25},
    'H': {'CAT': 0.41, 'CAC': 0.59}, 'I': {'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16},
    'L': {'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.40},
    'K': {'AAA': 0.42, 'AAG': 0.58}, 'M': {'ATG': 1.0}, 'F': {'TTT': 0.45, 'TTC': 0.55},
    'P': {'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.12},
    'S': {'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'AGT': 0.15, 'AGC': 0.24},
    'T': {'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12},
    'W': {'TGG': 1.0}, 'Y': {'TAT': 0.43, 'TAC': 0.57},
    'V': {'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47},
    '*': {'TAA': 0.30, 'TAG': 0.20, 'TGA': 0.50}
}
CODON_TO_AA = {codon: aa for aa, codons in CODON_USAGE_BIAS.items() for codon in codons}


# ────────────────────────────────────────────────────────
# 工具: 支持文件/目录的FASTA加载
# ────────────────────────────────────────────────────────
def iter_fasta(path):
    """兼容文件和目录的FASTA迭代器"""
    if os.path.isdir(path):
        for f in sorted(Path(path).glob("*.fasta")):
            for rec in SeqIO.parse(f, "fasta"):
                yield rec
    else:
        for rec in SeqIO.parse(path, "fasta"):
            yield rec

def load_genome_dict(path):
    """加载FASTA为 {acc: seq}"""
    genomes = {}
    for rec in iter_fasta(path):
        acc = rec.id.split()[0]
        if len(rec.seq) >= 200:
            genomes[acc] = str(rec.seq)
    return genomes

# ────────────────────────────────────────────────────────
# 双名法物种分类信息解析器
# ────────────────────────────────────────────────────────
def parse_fasta_species_categories(fasta_path):
    categories = {}
    for rec in iter_fasta(fasta_path):
        acc = rec.id.split()[0]
        desc = rec.description
        if desc.startswith(acc):
            desc = desc[len(acc):].strip()
        
        words = desc.split()
        if len(words) >= 2:
            species = f"{words[0]} {words[1]}"
        elif len(words) == 1:
            species = words[0]
        else:
            species = "Unknown Virus"
        categories[acc] = species
    return categories


# ────────────────────────────────────────────────────────
# 突变与反翻译底层引擎
# ────────────────────────────────────────────────────────
def mutate_viral_genome_biologically(seq, seq_id, total_rate, rng):
    if total_rate <= 0:
        return seq
    snv_rate = total_rate * 0.985
    ins_rate = total_rate * 0.0075
    del_rate = total_rate * 0.0075

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode='w') as infile, \
         tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode='r') as outfile:
        infile_path = infile.name
        outfile_path = outfile.name
        infile.write(f">{seq_id}\n{seq}\n")
        infile.close()
        outfile.close()

    try:
        cmd = [
            "mutation-simulator", "-i", infile_path, "-o", outfile_path,
            "--snv", f"{snv_rate:.6f}", "--ins", f"{ins_rate:.6f}", "--del", f"{del_rate:.6f}",
            "--dup", "0.0", "--inv", "0.0", "--tra", "0.0"
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        records = list(SeqIO.parse(outfile_path, "fasta"))
        mutated_seq = str(records[0].seq) if records else seq
    except (subprocess.CalledProcessError, FileNotFoundError):
        mutated_seq = fallback_mutate(seq, total_rate, rng)
    finally:
        if os.path.exists(infile_path): os.remove(infile_path)
        if os.path.exists(outfile_path): os.remove(outfile_path)
    return mutated_seq


def fallback_mutate(seq, rate, rng):
    s = list(seq)
    n = int(len(s) * rate)
    pos = rng.sample(range(len(s)), min(n, len(s)))
    for p in pos:
        orig = s[p]
        s[p] = rng.choice([b for b in 'ACGT' if b != orig])
    return ''.join(s)


def mutate_noncoding(seq, target_similarity, rng):
    return fallback_mutate(seq, 1.0 - target_similarity, rng)


def apply_safe_mutations(nt_seq, mutation_rate, rng):
    nt_list = list(nt_seq)
    seq_len = len(nt_seq)
    stop_codons = {'TAA', 'TAG', 'TGA'}
    n_mutations = int(seq_len * mutation_rate)
    positions = rng.sample(range(seq_len), min(n_mutations, seq_len))
    for p in positions:
        codon_idx = p // 3
        start = codon_idx * 3
        if start + 3 > seq_len: continue
        current_codon = ''.join(nt_list[start:start+3])
        p_in_codon = p % 3
        valid_bases = []
        for b in ['A', 'C', 'G', 'T']:
            if b == nt_list[p]: continue
            temp_codon = list(current_codon)
            temp_codon[p_in_codon] = b
            temp_codon_str = ''.join(temp_codon)
            if current_codon not in stop_codons:
                if temp_codon_str not in stop_codons: valid_bases.append(b)
            else:
                if temp_codon_str in stop_codons: valid_bases.append(b)
        if valid_bases:
            nt_list[p] = rng.choice(valid_bases)
    return ''.join(nt_list)


def backtranslate_to_similarity(orig_coding_seq, target_similarity, rng):
    codons = [orig_coding_seq[i:i+3] for i in range(0, len(orig_coding_seq) - 2, 3)]
    S_base = 0.68
    if target_similarity > S_base:
        p_keep = (target_similarity - S_base) / (1.0 - S_base)
        extra_mutation_rate = 0.0
    else:
        p_keep = 0.0
        extra_mutation_rate = (S_base - target_similarity) / S_base

    hybrid_codons = []
    for codon in codons:
        if codon not in CODON_TO_AA:
            hybrid_codons.append(codon)
            continue
        aa = CODON_TO_AA[codon]
        if rng.random() < p_keep:
            hybrid_codons.append(codon)
        else:
            choices = list(CODON_USAGE_BIAS[aa].keys())
            weights = list(CODON_USAGE_BIAS[aa].values())
            hybrid_codons.append(rng.choices(choices, weights=weights, k=1)[0])
    nt_seq = "".join(hybrid_codons)
    if extra_mutation_rate > 0.0:
        nt_seq = apply_safe_mutations(nt_seq, extra_mutation_rate, rng)
    return nt_seq


def process_genome_viromock(seq, target_similarity, rng, gf=None):
    if gf is None or not HAS_PYRODIGAL:
        # 降级: 整序列作为大CDS反向翻译 (nt~73%, aa~100%, 足够区分信号依赖)
        return backtranslate_to_similarity(seq, target_similarity, rng)
    genes = gf.find_genes(seq.encode('ascii'))
    sorted_genes = sorted(genes, key=lambda g: g.begin)
    reconstructed = []
    last_idx = 0
    for gene in sorted_genes:
        intergenic = seq[last_idx:gene.begin]
        if intergenic:
            reconstructed.append(mutate_noncoding(intergenic, target_similarity, rng))
        gene_seq = seq[gene.begin:gene.end]
        if gene.strand == -1:
            gene_seq = str(Seq(gene_seq).reverse_complement())
        coding_nt = backtranslate_to_similarity(gene_seq, target_similarity, rng)
        if gene.strand == -1:
            coding_nt = str(Seq(coding_nt).reverse_complement())
        reconstructed.append(coding_nt)
        last_idx = gene.end
    if last_idx < len(seq):
        reconstructed.append(mutate_noncoding(seq[last_idx:], target_similarity, rng))
    return "".join(reconstructed)


# ────────────────────────────────────────────────────────
# 裁剪并生成标准竖线标头ID
# ────────────────────────────────────────────────────────
def extract_fragments_pipe_format(seq, cov_pct, n_rep, scheme_type, mut_or_sim_val, source_acc, rng):
    fragments = []
    L = len(seq)
    frag_len = max(200, int(L * cov_pct / 100))
    frag_len = min(frag_len, L)
    
    for r in range(1, n_rep + 1):
        start = rng.randint(0, max(1, L - frag_len)) if L > frag_len else 0
        frag_seq = seq[start:start + frag_len]
        
        if scheme_type in ["A", "B"]:
            seq_id = f"positive|scheme_{scheme_type}|mut{mut_or_sim_val}|cov{cov_pct}|rep{r}|source={source_acc}"
        else:
            seq_id = f"positive|scheme_C|sim{mut_or_sim_val}|cov100|rep{r}|source={source_acc}"
            
        fragments.append({
            "seq_id": seq_id,
            "seq": frag_seq,
            "source": source_acc,
            "coverage": cov_pct,
            "rep": r
        })
    return fragments


# ────────────────────────────────────────────────────────
# 主程序
# ────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="MASTER版评估数据集一键构建流水线 (含本地序列注入功能)")
    parser.add_argument("--known-fasta", required=True, help="已知病毒FASTA (60条)")
    parser.add_argument("--novel-fasta", required=True, help="标准候选新病毒FASTA (45条)")
    parser.add_argument("--extra-novel-fasta", default=None,
                        help="额外追加的本地新病毒 FASTA 文件 (可选，合并到B组和C组)")
    parser.add_argument("--virus-meta", default=None,
                        help="selected_viruses.tsv (含 accession/species 列, 用于Category映射)")
    parser.add_argument("--conserved-fasta", required=True, 
                        help="保守结构域 FASTA 文件，或者包含多个 trap_*.fasta 的目录")
    parser.add_argument("--eve-fasta", required=True, help="内源性病毒元件 FASTA")
    parser.add_argument("--host-fasta", required=True, help="宿主基因组 FASTA")
    parser.add_argument("--outdir", required=True, help="输出目录")
    
    # 终极优化科学参数
    parser.add_argument("--known-mutations", type=int, nargs="+", default=[0, 5, 10], help="方案A(已知)的突变梯度")
    parser.add_argument("--similarity", type=int, default=70, help="方案C(VIROMOCK)目标相似度")
    parser.add_argument("--coverage-levels", type=int, nargs="+", default=[100, 90, 70, 50, 30, 10], help="A/B方案覆盖率裁剪梯度")
    parser.add_argument("--n-per-coverage", type=int, default=5, help="裁剪重复次数与方案C独立回译次数")
    parser.add_argument("--n-neg", type=int, default=1000, help="每类负样本(D/E/H)的数量")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    gf = pyrodigal.GeneFinder(meta=True) if HAS_PYRODIGAL else None

    # 1. 从 selected_viruses.tsv 加载 accession→species 映射
    acc2species = {}
    if args.virus_meta and os.path.exists(args.virus_meta):
        meta = pd.read_csv(args.virus_meta, sep='\t')
        acc2species = dict(zip(meta['accession'], meta['species']))
        print(f"[Meta] Loaded {len(acc2species)} accessions from {args.virus_meta}")
    # 补充：FASTA头解析 (仅用于不在meta中的新病毒)
    def get_species(acc, fasta_cats=None):
        if acc in acc2species:
            return acc2species[acc]
        if fasta_cats and acc in fasta_cats and fasta_cats[acc] != 'Unknown Virus':
            return fasta_cats[acc]
        return acc  # 最后用accession本身

    print("[Parse] Assigning species from metadata...")
    known_cats_from_fasta = parse_fasta_species_categories(args.known_fasta)
    novel_cats_from_fasta = parse_fasta_species_categories(args.novel_fasta)

    # 2. 隔离加载已知/新病毒基因组
    known_genomes = load_genome_dict(args.known_fasta)
    novel_genomes = load_genome_dict(args.novel_fasta)

    # 3. 智能合并追加的本地新病毒
    if args.extra_novel_fasta and os.path.exists(args.extra_novel_fasta):
        print(f"[Parse] Merging extra local novel genomes from: {args.extra_novel_fasta}")
        extra_cats = parse_fasta_species_categories(args.extra_novel_fasta)
        novel_cats_from_fasta.update(extra_cats)
        extra_genomes = load_genome_dict(args.extra_novel_fasta)
        novel_genomes.update(extra_genomes)
        print(f"[Parse] Successfully added {len(extra_genomes)} local genomes to Novel Pool.")

    # 构建全局分配函数
    def assign_category(acc):
        """为每个accession分配Category (优先metadata, 其次FASTA头解析, 最后用accession)"""
        return get_species(acc, {**known_cats_from_fasta, **novel_cats_from_fasta})
    all_genomes = {**known_genomes, **novel_genomes}
    
    print(f"[Init] Active Genomes: {len(known_genomes)} Known | {len(novel_genomes)} Novel (standard + extra)")

    # 汇总数据
    master_labels = []
    all_eval_records = []

    # ==========================================
    # 方案 A: 已知病毒评估 (仅使用已知病毒)
    # ==========================================
    print(f"\n[Scheme A] Simulating Known Evaluation (muts: {args.known_mutations})...")
    for acc, seq in known_genomes.items():
        species_cat = assign_category(acc)
        for mut in args.known_mutations:
            mut_seq = mutate_viral_genome_biologically(seq, acc, mut/100.0, rng)
            for cov in args.coverage_levels:
                frags = extract_fragments_pipe_format(mut_seq, cov, args.n_per_coverage, "A", mut, acc, rng)
                for f in frags:
                    all_eval_records.append(SeqRecord(Seq(f["seq"]), id=f["seq_id"], description=""))
                    master_labels.append({
                        "seq_id": f["seq_id"],
                        "label": "positive",
                        "type": "scheme_A",
                        "Category": species_cat
                    })

    # ==========================================
    # 方案 B: 真实新病毒评估 (标准新病毒 + 追加的本地序列，严格 0%)
    # ==========================================
    print(f"\n[Scheme B] Simulating Novel Evaluation (0% Wild-type only)...")
    for acc, seq in novel_genomes.items():
        species_cat = assign_category(acc)
        for cov in args.coverage_levels:
            frags = extract_fragments_pipe_format(seq, cov, args.n_per_coverage, "B", 0, acc, rng)
            for f in frags:
                all_eval_records.append(SeqRecord(Seq(f["seq"]), id=f["seq_id"], description=""))
                master_labels.append({
                    "seq_id": f["seq_id"],
                    "label": "positive",
                    "type": "scheme_B",
                    "Category": species_cat
                })

    # ==========================================
    # 方案 C: VIROMOCK 评估 (全部 105+ 病毒，100% 全长)
    # ==========================================
    print(f"\n[Scheme C] Simulating VIROMOCK synonymous evaluation (100% Full-Length only)...")
    sim_float = args.similarity / 100.0
    for acc, seq in all_genomes.items():
        species_cat = assign_category(acc)
        for r in range(1, args.n_per_coverage + 1):
            mut_seq = process_genome_viromock(seq, sim_float, rng, gf)
            frags = extract_fragments_pipe_format(mut_seq, 100, 1, "C", args.similarity, acc, rng)
            f = frags[0]
            corrected_id = f"positive|scheme_C|sim{args.similarity}|cov100|rep{r}|source={acc}"
            
            all_eval_records.append(SeqRecord(Seq(f["seq"]), id=corrected_id, description=""))
            master_labels.append({
                "seq_id": corrected_id,
                "label": "positive",
                "type": "scheme_C",
                "Category": species_cat
            })

    # ==========================================
    # 阴性对照组提取 (D, E, H 各 n_neg 条)
    # ==========================================
    print(f"\n[Negatives] Extracting background controls (D/E/H each {args.n_neg})...")
    all_frag_lens = [len(r.seq) for r in all_eval_records[:2000]] if all_eval_records else [500]

    # D: 保守结构域陷阱
    d_seqs = []
    path_trap = Path(args.conserved_fasta)
    if path_trap.is_dir():
        print(f"  [Loader] Loading conserved traps from directory: {path_trap}")
        for f in sorted(path_trap.glob("*.fasta")):
            for rec in SeqIO.parse(f, "fasta"):
                if len(rec.seq) >= 200: d_seqs.append(str(rec.seq))
    elif path_trap.is_file():
        print(f"  [Loader] Loading conserved traps from single file: {path_trap}")
        for rec in SeqIO.parse(path_trap, "fasta"):
            if len(rec.seq) >= 200: d_seqs.append(str(rec.seq))
    else:
        print(f"  [WARNING] --conserved-fasta target '{path_trap}' not found. Using synthetic fallback.")
        d_seqs = ["ACGT" * 500]

    # 生成 D 组
    for idx in range(args.n_neg):
        ref_seq = rng.choice(d_seqs)
        flen = min(rng.choice(all_frag_lens), len(ref_seq))
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(ref_seq) - flen))
        frag = ref_seq[start:start+flen]
        
        seq_id = f"negative|conserved_trap|idx{idx:04d}|source=conserved_trap"
        all_eval_records.append(SeqRecord(Seq(frag), id=seq_id, description=""))
        master_labels.append({
            "seq_id": seq_id,
            "label": "negative",
            "type": "conserved_trap",
            "Category": "conserved_trap"
        })

    # E: 内源性病毒元件
    e_seqs = [str(r.seq) for r in SeqIO.parse(args.eve_fasta, "fasta") if len(r.seq) >= 200]
    for idx in range(args.n_neg):
        ref_seq = rng.choice(e_seqs)
        flen = min(rng.choice(all_frag_lens), len(ref_seq))
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(ref_seq) - flen))
        frag = ref_seq[start:start+flen]
        
        seq_id = f"negative|EVE|idx{idx:04d}|source=EVE_transposon"
        all_eval_records.append(SeqRecord(Seq(frag), id=seq_id, description=""))
        master_labels.append({
            "seq_id": seq_id,
            "label": "negative",
            "type": "EVE",
            "Category": "EVE_transposon"
        })

    # H: 宿主基因组随机噪音
    h_seqs = [str(r.seq) for r in SeqIO.parse(args.host_fasta, "fasta") if len(r.seq) >= 5000]
    for idx in range(args.n_neg):
        ref_seq = rng.choice(h_seqs)
        flen = min(rng.choice(all_frag_lens), len(ref_seq) - 10)
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(ref_seq) - flen))
        frag = ref_seq[start:start+flen]
        
        seq_id = f"negative|host_random|idx{idx:04d}|source=host_genomic_background"
        all_eval_records.append(SeqRecord(Seq(frag), id=seq_id, description=""))
        master_labels.append({
            "seq_id": seq_id,
            "label": "negative",
            "type": "host_random",
            "Category": "host_genomic_background"
        })

    # ==========================================
    # 写出终极评估数据集与主标签TSV
    # ==========================================
    fasta_out = os.path.join(args.outdir, "evaluation_sequences.fasta")
    labels_out = os.path.join(args.outdir, "sequence_labels_category.tsv")

    SeqIO.write(all_eval_records, fasta_out, "fasta")
    pd.DataFrame(master_labels).to_csv(labels_out, sep='\t', index=False)

    print(f"\n{'='*60}\n   MASTER 评估数据集一键构建成功！")
    print(f"  输出路径: {args.outdir}")
    print(f"  合并 FASTA: {fasta_out}")
    print(f"  主分类标签 TSV:  {labels_out}")
    print(f"  总序列数: {len(all_eval_records)} 条")
    print(f"    - 阳性 A (已知):   {len([l for l in master_labels if l['type'] == 'scheme_A'])} 条")
    print(f"    - 阳性 B (新种):   {len([l for l in master_labels if l['type'] == 'scheme_B'])} 条")
    print(f"    - 阳性 C (同义):   {len([l for l in master_labels if l['type'] == 'scheme_C'])} 条")
    print(f"    - 阴性背景 (D/E/H): {args.n_neg * 3} 条")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()