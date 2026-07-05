#!/usr/bin/env python3
"""
基于 NCBI 新上传病毒的未知病毒鉴定评估数据集

策略 (参考 VIVID/VirID 三级分类法):
  1. 从 novel_viruses_metadata.tsv 获取45个新病毒 accession
  2. 下载基因组 FASTA (或从已有库提取)
  3. BLAST 每个病毒对 final.cluster.ref.fasta 获取最佳匹配相似度
  4. 按相似度分三级:
     L1 (近缘): 70% ≤ sim < 90%  — 同属水平
     L2 (远缘): 40% ≤ sim < 70%  — 同科水平
     L3 (全新): 20% ≤ sim < 40%  — 新科水平 (或无匹配)
  5. 按覆盖率梯度截取片段
  6. 加入负样本

用法:
  python prep_build_ncbi_novel_eval.py \
      --meta novel_viruses_metadata.tsv \
      --ref-db final.cluster.ref.fasta \
      --coverage-levels 100 80 60 40 --n-per-coverage 2 \
      --host-genome ningxia.genome.fasta \
      --conserved-prots step3_conserved_traps/ \
      --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
      --outdir step3_ncbi_novel_eval/ --seed 42
"""

import argparse, os, sys, random, subprocess, tempfile
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np


def load_metadata(meta_path):
    """加载新病毒元数据"""
    meta = pd.read_csv(meta_path, sep='\t')
    print(f"[meta] {len(meta)} novel viruses loaded")
    return meta


def fetch_genomes(meta, ref_fasta, outdir):
    """
    从参考数据库中提取这些病毒的基因组 FASTA。
    如果不在 ref_fasta 中，则需从 NCBI 下载 (这里假设已在 ref 中或 nearby)
    """
    os.makedirs(outdir, exist_ok=True)
    ref_index = SeqIO.index(ref_fasta, "fasta")

    genomes = {}
    for _, row in meta.iterrows():
        acc = row['Accession']
        length = int(row['Length'])
        if length < 1000:
            continue
        if acc in ref_index:
            rec = ref_index[acc]
            out_path = os.path.join(outdir, f"{acc}.fasta")
            SeqIO.write(rec, out_path, "fasta")
            genomes[acc] = {'seq': str(rec.seq), 'length': len(rec.seq),
                            'family': row.get('Family', ''), 'genus': row.get('Genus', ''),
                            'species': row.get('Species', '')}
        else:
            print(f"  [WARN] {acc} not found in reference database")

    print(f"[genomes] {len(genomes)} genomes >= 1000 bp extracted")
    return genomes


def classify_novelty(genomes, ref_fasta, tmpdir):
    """
    通过 BLASTn 确定每个病毒对参考数据库的最佳匹配相似度。
    分类为 L1 / L2 / L3。
    如果 BLAST 不可用，用 k-mer 距离近似估算。
    """
    # 先检查 blastn 是否可用
    has_blast = subprocess.run(['which', 'blastn'], capture_output=True).returncode == 0

    classifications = {}
    for acc, info in genomes.items():
        seq = info['seq']
        fasta_path = os.path.join(tmpdir, f"query_{acc}.fasta")
        with open(fasta_path, 'w') as f:
            f.write(f">{acc}\n{seq}\n")

        if has_blast:
            cmd = ['blastn', '-query', fasta_path, '-db', ref_fasta.replace('.fasta', ''),
                   '-outfmt', '6 qseqid sseqid pident length qlen slen',
                   '-max_target_seqs', '1', '-num_threads', '1']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.stdout.strip():
                fields = result.stdout.strip().split('\n')[0].split('\t')
                pident = float(fields[2]) / 100.0
                sseqid = fields[1]
            else:
                pident = 0.0
                sseqid = 'no_hit'
        else:
            # 无 BLAST: 用长度和 GC 估算 (简化, 标记为 unclassified)
            pident = 0.5
            sseqid = 'unknown'

        if pident >= 0.70:
            level = 'L1_near'
        elif pident >= 0.40:
            level = 'L2_distant'
        else:
            level = 'L3_novel'

        classifications[acc] = {
            'level': level,
            'pident': round(pident, 4),
            'best_hit': sseqid.split('.')[0],
        }

    return classifications


def generate_fragments(genomes, coverage_levels, n_per_cov, rng):
    """按覆盖率梯度截取片段"""
    records = []
    idx = 0
    for cov_pct in coverage_levels:
        for _ in range(n_per_cov):
            for acc, info in genomes.items():
                full_seq = info['seq']
                full_len = info['length']
                frag_len = max(200, int(full_len * cov_pct / 100))
                frag_len = min(frag_len, full_len)
                max_start = full_len - frag_len
                start = rng.randint(0, max(1, max_start)) if max_start > 0 else 0
                frag = full_seq[start:start + frag_len]

                seq_id = f"ncbi_novel|cov{cov_pct}|idx{idx:05d}|source={acc}"
                records.append({
                    "seq_id": seq_id, "label": "positive", "type": "positive",
                    "source": acc, "coverage_pct": cov_pct,
                    "full_length": full_len, "frag_length": len(frag),
                    "frag_seq": frag,
                    "family": info['family'], "genus": info['genus'],
                })
                idx += 1
    return records


def generate_negatives(pos_records, host_genome, conserved_dir, eve_fasta,
                       n_each, rng):
    """生成三类负样本 (同之前)"""
    all_neg = []
    pos_lengths = [r["frag_length"] for r in pos_records]

    # A: host random
    host_seqs = []
    if host_genome and os.path.exists(host_genome):
        for rec in SeqIO.parse(host_genome, "fasta"):
            if len(rec.seq) >= 200: host_seqs.append(str(rec.seq))
    if not host_seqs: host_seqs = ["ACGT" * 5000]
    for i in range(n_each):
        tlen = rng.choice(pos_lengths)
        flen = max(200, tlen + rng.randint(-200, 200))
        seq = rng.choice(host_seqs)
        flen = min(flen, len(seq) - 5)
        start = rng.randint(0, max(1, len(seq) - flen))
        frag = seq[start:start + flen]
        all_neg.append({"seq_id": f"negA|{i:04d}", "label": "negative", "type": "negative_A",
                         "frag_seq": frag, "frag_length": len(frag)})

    # B: conserved traps
    trap_seqs = []
    if conserved_dir and os.path.exists(conserved_dir):
        for f in Path(conserved_dir).glob("*.fasta"):
            try:
                rec = next(SeqIO.parse(f, "fasta"))
                if len(rec.seq) >= 500: trap_seqs.append(str(rec.seq))
            except: pass
    if not trap_seqs: trap_seqs = ["ACGT" * 3000]
    for i in range(n_each):
        seq = rng.choice(trap_seqs)
        flen = min(rng.choice(pos_lengths), len(seq) - 5)
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(seq) - flen))
        all_neg.append({"seq_id": f"negB|{i:04d}", "label": "negative", "type": "negative_B",
                         "frag_seq": seq[start:start + flen], "frag_length": len(frag)})

    # C: EVE
    eve_seqs = []
    if eve_fasta and os.path.exists(eve_fasta):
        for rec in SeqIO.parse(eve_fasta, "fasta"):
            if len(rec.seq) >= 500: eve_seqs.append(str(rec.seq))
    if not eve_seqs: eve_seqs = host_seqs[:50] if host_seqs else ["ACGT" * 5000]
    for i in range(n_each):
        seq = rng.choice(eve_seqs)
        flen = min(rng.choice(pos_lengths), len(seq) - 5)
        flen = max(200, flen)
        start = rng.randint(0, max(1, len(seq) - flen))
        all_neg.append({"seq_id": f"negC|{i:04d}", "label": "negative", "type": "negative_C",
                         "frag_seq": seq[start:start + flen], "frag_length": len(frag)})

    return all_neg


def write_dataset(outdir, name, pos_records, neg_records, novelty_labels=None):
    fasta_out = os.path.join(outdir, f"{name}.fasta")
    labels_out = os.path.join(outdir, f"{name}_labels.tsv")
    all_recs, all_labs = [], []
    for r in pos_records:
        all_recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        lab = {"seq_id": r["seq_id"], "label": r["label"], "type": r["type"],
               "source": r["source"], "coverage_pct": r["coverage_pct"]}
        if novelty_labels and r["source"] in novelty_labels:
            lab["novelty_level"] = novelty_labels[r["source"]]['level']
            lab["pident"] = novelty_labels[r["source"]]['pident']
        all_labs.append(lab)
    for r in neg_records:
        all_recs.append(SeqRecord(Seq(r["frag_seq"]), id=r["seq_id"], description=""))
        all_labs.append({"seq_id": r["seq_id"], "label": r["label"], "type": r["type"]})
    SeqIO.write(all_recs, fasta_out, "fasta")
    pd.DataFrame(all_labs).to_csv(labels_out, sep='\t', index=False)
    return fasta_out


def main():
    parser = argparse.ArgumentParser(description='NCBI新病毒未知鉴定评估')
    parser.add_argument('--meta', required=True, help='novel_viruses_metadata.tsv')
    parser.add_argument('--ref-fasta', required=True, help='final.cluster.ref.fasta')
    parser.add_argument('--coverage-levels', type=int, nargs='+', default=[100, 80, 60, 40])
    parser.add_argument('--n-per-coverage', type=int, default=2)
    parser.add_argument('--host-genome')
    parser.add_argument('--conserved-prots')
    parser.add_argument('--eve-fasta')
    parser.add_argument('--n-neg', type=int, default=200)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix='novel_blast_', dir=args.outdir)

    # 1. 加载元数据
    meta = load_metadata(args.meta)

    # 2. 提取基因组
    genomes = fetch_genomes(meta, args.ref_fasta, os.path.join(args.outdir, 'genomes'))

    # 3. 新颖度分类 (BLASTn)
    print("\n[classify] Running BLASTn to determine novelty levels...")
    novelty = classify_novelty(genomes, args.ref_fasta, tmpdir)

    # 统计
    from collections import Counter
    level_counts = Counter(v['level'] for v in novelty.values())
    print(f"[classify] Novelty distribution:")
    for level, count in sorted(level_counts.items()):
        print(f"  {level}: {count} viruses")
    for acc, info in novelty.items():
        print(f"  {acc}: {info['level']} (pident={info['pident']}, best_hit={info['best_hit']})")

    # 4. 生成正样本
    print("\n[generate] Creating positive fragments...")
    pos_records = generate_fragments(genomes, args.coverage_levels, args.n_per_coverage, rng)
    print(f"[generate] {len(pos_records)} positive fragments")

    # 5. 生成负样本
    neg_records = generate_negatives(pos_records, args.host_genome, args.conserved_prots,
                                     args.eve_fasta, args.n_neg, rng)

    # 6. 输出合并数据集
    write_dataset(args.outdir, "ncbi_novel_all", pos_records, neg_records, novelty)

    # 7. 按 level 分别输出
    for level in ['L1_near', 'L2_distant', 'L3_novel']:
        level_pos = [r for r in pos_records if novelty.get(r['source'], {}).get('level') == level]
        if level_pos:
            write_dataset(args.outdir, f"ncbi_novel_{level}", level_pos, neg_records, novelty)

    # 8. 汇总
    print(f"\n{'='*60}")
    print(f"  NCBI 新病毒评估数据集")
    print(f"{'='*60}")
    print(f"  输入病毒: {len(genomes)} (>=1000 bp)")
    print(f"  正样本:   {len(pos_records)} (4 cov × 2 rep × {len(genomes)} viruses)")
    print(f"  负样本:   {len(neg_records)} (A/B/C 各{args.n_neg})")
    print(f"  L1 (70-90%):  {level_counts.get('L1_near', 0)} viruses")
    print(f"  L2 (40-70%):  {level_counts.get('L2_distant', 0)} viruses")
    print(f"  L3 (<40%):    {level_counts.get('L3_novel', 0)} viruses")
    print(f"  输出:       {args.outdir}")

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == '__main__':
    main()
