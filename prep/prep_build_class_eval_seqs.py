#!/usr/bin/env python3
"""
构建病毒分类方法评估数据集（v2 — 覆盖度梯度 + 防泄漏）

从同一批50个评估用病毒（或者单独指定的ICTV参考病毒），按覆盖度梯度截取测试序列。
实现"防信息泄漏"——每个测试序列的同种序列从参考数据库中排除。

用法:
  # 模式A：使用同一批50个评估用病毒
  python prep_build_class_eval_seqs.py \
      --virus-dir eval_viruses_50/ \
      --ref-info final.cluster.ref_info.tsv \
      --ref-fasta final.cluster.ref.fasta \
      --coverage-levels 100 80 60 40 20 \
      --n-per-coverage 2 \
      --outdir eval_classification/ --seed 42

  # 模式B：从VMR中选择独立测试病毒
  python prep_build_class_eval_seqs.py \
      --vmr VMR_MSL41.v2.tsv \
      --ref-info final.cluster.ref_info.tsv \
      --ref-fasta final.cluster.ref.fasta \
      --n-total 300 \
      --outdir eval_classification/ --seed 42

输出:
  {outdir}/
    ├── test_sequences/             # 测试序列（按覆盖度梯度截取）
    ├── test_metadata.tsv           # 真值（科/属/种/覆盖度/长度）
    ├── db_sequences.fasta          # 去泄漏参考数据库（单个FASTA）
    └── leakage_report.tsv          # 信息泄漏检查报告
"""

import argparse, os, sys, random
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd


def load_ref_info(path):
    """加载 ref_info TSV → {accession: {field: value}}"""
    import csv
    info = {}
    with open(path, "r") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc = row.get("Accession", "").strip()
            if acc:
                info[acc] = {
                    "species": row.get("VMR_Species", "").strip(),
                    "genus": row.get("VMR_Genus", "").strip(),
                    "family": row.get("VMR_Family", "").strip(),
                }
    print(f"[ref_info] Loaded {len(info)} records")
    return info


def load_virus_sequences(virus_input, ref_info):
    """
    加载病毒序列。
    virus_input 可以是目录(*.fasta)或单个FASTA文件。
    返回 [(accession, seq, species, genus, family, length), ...]
    """
    seqs = []
    if os.path.isdir(virus_input):
        files = sorted(Path(virus_input).glob("*.fasta"))
    else:
        files = [Path(virus_input)]

    for f in files:
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            info = ref_info.get(acc, {})
            species = info.get("species", "")
            genus = info.get("genus", "")
            family = info.get("family", "")
            if not species or not genus or not family:
                continue
            seqs.append((acc, str(rec.seq), species, genus, family, len(rec.seq)))

    print(f"[virus] Loaded {len(seqs)} sequences with complete ICTV classification")
    return seqs


def generate_test_sequences_by_coverage(virus_seqs, coverage_levels, n_per_cov, rng):
    """
    按覆盖度梯度从病毒全长基因组截取测试序列。
    返回 test_records + cut_seqs 列表。
    """
    test_records = []
    cut_seqs = []
    idx = 0

    for cov_pct in coverage_levels:
        for _ in range(n_per_cov):
            for acc, full_seq, species, genus, family, full_len in virus_seqs:
                frag_len = max(500, int(full_len * cov_pct / 100))
                if frag_len > full_len:
                    frag_len = full_len

                max_start = full_len - frag_len
                start = rng.randint(0, max(1, max_start)) if max_start > 0 else 0
                frag = full_seq[start:start + frag_len]

                seq_id = f"test|cov{cov_pct}|{idx:04d}|src={acc}"
                cut_seqs.append((seq_id, frag))
                test_records.append({
                    "seq_id": seq_id,
                    "source_accession": acc,
                    "species": species, "genus": genus, "family": family,
                    "coverage_pct": cov_pct,
                    "full_length": full_len,
                    "frag_length": len(frag),
                })
                idx += 1

    print(f"[test] Generated {len(test_records)} test sequences "
          f"({len(virus_seqs)} viruses × {len(coverage_levels)} coverage × {n_per_cov})")
    return test_records, cut_seqs


def build_leakage_free_db(virus_input, ref_info, test_species, test_accessions):
    """
    构建去泄漏参考数据库：排除测试序列的 accession 和同 species。
    写入单个 FASTA 文件。
    """
    if os.path.isdir(virus_input):
        files = sorted(Path(virus_input).glob("*.fasta"))
    else:
        files = [Path(virus_input)]

    excluded_acc = 0
    excluded_sp = 0
    kept = 0

    db_records = []
    for f in files:
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            if acc in test_accessions:
                excluded_acc += 1
                continue
            info = ref_info.get(acc, {})
            sp = info.get("species", "")
            if sp and sp in test_species:
                excluded_sp += 1
                continue
            db_records.append(rec)
            kept += 1

    print(f"[db] Kept {kept} sequences (excluded {excluded_acc} exact + {excluded_sp} same-species)")

    if excluded_acc + excluded_sp == 0:
        print("[db] WARNING: No sequences excluded. Check species name matching between "
              "test metadata and ref_info.")

    return db_records


def main():
    parser = argparse.ArgumentParser(description="构建病毒分类评估数据集 (v2)")
    parser.add_argument("--virus-dir", help="评估用病毒目录（FASTA文件）")
    parser.add_argument("--vmr", help="ICTV VMR TSV (模式B)")
    parser.add_argument("--ref-info", required=True, help="final.cluster.ref_info.tsv")
    parser.add_argument("--ref-fasta", required=True, help="final.cluster.ref.fasta")
    parser.add_argument("--n-total", type=int, default=300,
                        help="模式B的目标测试序列数（模式A由coverage参数控制）")
    parser.add_argument("--coverage-levels", type=int, nargs="+",
                        default=[100, 80, 60, 40, 20])
    parser.add_argument("--n-per-coverage", type=int, default=2,
                        help="每个病毒每种覆盖度的截取次数")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    ref_info = load_ref_info(args.ref_info)

    if args.virus_dir:
        # 模式A：使用评估用病毒目录
        virus_seqs = load_virus_sequences(args.virus_dir, ref_info)
    else:
        # 模式B：从 ref_fasta + VMR 中筛选
        virus_seqs = load_virus_sequences(args.ref_fasta, ref_info)
        if args.vmr and os.path.exists(args.vmr):
            # 进一步按 VMR 的 ICTV 分类筛选
            vmr = pd.read_csv(args.vmr, sep="\t", low_memory=False)
            vmr_species = set(str(s).strip() for s in vmr.get("Species", []))
            virus_seqs = [v for v in virus_seqs if v[2] in vmr_species]
            print(f"[vmr] Filtered to {len(virus_seqs)} sequences in VMR")

    if len(virus_seqs) < 5:
        print("[ERROR] Too few virus sequences with complete classification!")
        sys.exit(1)

    # 生成测试序列
    test_records, cut_seqs = generate_test_sequences_by_coverage(
        virus_seqs, args.coverage_levels, args.n_per_coverage, rng)

    # 写出测试序列
    test_dir = os.path.join(args.outdir, "test_sequences")
    os.makedirs(test_dir, exist_ok=True)
    for seq_id, seq_str in cut_seqs:
        out_file = os.path.join(test_dir, f"{seq_id.replace('|', '_')}.fasta")
        with open(out_file, "w") as f:
            f.write(f">{seq_id}\n{seq_str}\n")

    # 写出测试元数据
    meta_df = pd.DataFrame(test_records)
    meta_df.to_csv(os.path.join(args.outdir, "test_metadata.tsv"), sep="\t", index=False)

    # 构建去泄漏数据库
    test_species = set(r["species"] for r in test_records)
    test_accessions = set(r["source_accession"] for r in test_records)
    db_records = build_leakage_free_db(
        args.virus_dir or args.ref_fasta, ref_info, test_species, test_accessions)

    db_fasta = os.path.join(args.outdir, "db_sequences.fasta")
    SeqIO.write(db_records, db_fasta, "fasta")

    # 统计
    cov_dist = Counter(r["coverage_pct"] for r in test_records)
    families = set(r["family"] for r in test_records)
    genera = set(r["genus"] for r in test_records)
    print(f"\n[DONE]")
    print(f"  Test sequences: {len(test_records)}")
    print(f"  Coverage: {dict(sorted(cov_dist.items()))}")
    print(f"  Families: {len(families)}  Genera: {len(genera)}")
    print(f"  DB sequences (leakage-free): {len(db_records)}")


if __name__ == "__main__":
    main()
