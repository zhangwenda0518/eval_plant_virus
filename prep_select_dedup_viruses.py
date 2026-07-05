#!/usr/bin/env python3
"""
去冗余评估用病毒选取脚本
特点：
 - 从 final.cluster.ref_info.tsv 按属(VMR_Genus)分组
 - 每属随机选取 3-5 个物种
 - 每个物种选一个代表性 accession（优先 longer genome）
 - 输出 FASTA 文件和选取记录

用法:
  python prep_select_dedup_viruses.py \
      --ref-info final.cluster.ref_info.tsv \
      --ref-fasta final.cluster.ref.old.fasta \
      --n-per-genus 3 5 \
      --outdir step1_dedup_viruses/ --seed 42
"""

import argparse, os, random, csv
from collections import defaultdict
from Bio import SeqIO
import pandas as pd


def load_ref_info(path):
    """加载 ref_info，筛选完整非节段病毒，按属分组"""
    by_genus = defaultdict(list)
    excluded = defaultdict(int)
    total = 0

    with open(path, "r") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            # 硬约束
            cat = str(row.get("Category", "")).strip().lower()
            if "segmented" in cat and "nonsegmented" not in cat:
                excluded["segmented"] += 1
                continue

            comp = str(row.get("Nuc_Completeness", "")).strip().lower()
            if comp != "complete":
                excluded["not_complete"] += 1
                continue

            genus = str(row.get("VMR_Genus", "")).strip()
            species = str(row.get("VMR_Species", "")).strip()
            acc = str(row.get("Accession", "")).strip()
            length_str = str(row.get("Length", "0")).strip()

            if not genus or not species or not acc:
                excluded["missing_meta"] += 1
                continue

            length = int(length_str) if length_str.isdigit() else 0
            host = str(row.get("Host", "")).strip()
            title = str(row.get("GenBank_Title", "")).strip().lower()
            # 排除非植物宿主（绿藻病毒等）
            non_plant_kw = ["chlorella", "phycodnaviridae", "chlorovirus",
                           "prasinovirus", "coccolithovirus", "prymnesiovirus",
                           "phaeovirus", "raphidovirus"]
            skip = any(kw in host.lower() or kw in title for kw in non_plant_kw)
            if skip:
                excluded["non_plant_host"] += 1
                continue

            by_genus[genus].append({
                "accession": acc, "species": species, "genus": genus,
                "length": length, "host": host,
                "family": str(row.get("VMR_Family", "")).strip(),
            })

    print(f"[ref_info] Total records: {total}")
    print(f"[ref_info] Excluded: {dict(excluded)}")
    print(f"[ref_info] Retained: {sum(len(v) for v in by_genus.values())} in {len(by_genus)} genera")
    return by_genus


def select_per_genus(by_genus, min_n, max_n, n_genera, rng):
    """每属选 min_n~max_n 个物种，n_genera=0 选全部属"""
    selected = []
    log = []

    genera_list = sorted(by_genus.keys())
    if n_genera > 0 and n_genera < len(genera_list):
        genera_list = rng.sample(genera_list, n_genera)
        print(f"[select] Limiting to {n_genera}/{len(by_genus)} genera")

    for genus in genera_list:
        members = by_genus[genus]
        # 按物种分组
        by_species = defaultdict(list)
        for m in members:
            by_species[m["species"]].append(m)

        species_list = list(by_species.keys())
        n_wanted = min(max_n, max(min_n, min(len(species_list), max_n)))
        if len(species_list) < min_n:
            print(f"  [SKIP] {genus}: only {len(species_list)} species (< {min_n})")
            continue

        chosen_species = rng.sample(species_list, n_wanted)
        for sp in chosen_species:
            # 选该物种中最长的 accession
            best = max(by_species[sp], key=lambda x: x["length"])
            selected.append(best)
            log.append(best)

        print(f"  {genus}: {len(species_list)} species → picked {n_wanted}")

    print(f"\n[select] Total: {len(selected)} viruses from {len(set(r['genus'] for r in log))} genera")
    return selected, log


def write_output(selected, log, fasta_path, outdir):
    """写出 FASTA 和选取记录"""
    os.makedirs(outdir, exist_ok=True)

    # 建立 FASTA 索引
    fasta_idx = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        acc = rec.id.split()[0]
        fasta_idx[acc] = rec
        fasta_idx[rec.id] = rec
    print(f"[fasta] Indexed {len(fasta_idx)//2} sequences")

    written = 0
    for item in selected:
        acc = item["accession"]
        rec = fasta_idx.get(acc)
        if rec is None:
            # 尝试模糊匹配
            for k, v in fasta_idx.items():
                if k.startswith(acc) or acc in k:
                    rec = v
                    break
        if rec is None:
            print(f"[WARNING] {acc} not found in FASTA")
            continue
        out_path = os.path.join(outdir, f"{acc}.fasta")
        SeqIO.write(rec, out_path, "fasta")
        written += 1

    df = pd.DataFrame(log)
    df.to_csv(os.path.join(outdir, "selected_viruses.tsv"), sep="\t", index=False)

    print(f"\n[DONE] Wrote {written} FASTA files to {outdir}")
    print(f"  Genera: {df['genus'].nunique()}")
    print(f"  Species: {df['species'].nunique()}")
    print(f"  Length range: {df['length'].min()} - {df['length'].max()} bp")


def main():
    parser = argparse.ArgumentParser(description="去冗余评估用病毒选取")
    parser.add_argument("--ref-info", required=True)
    parser.add_argument("--ref-fasta", required=True)
    parser.add_argument("--n-per-genus", type=int, nargs=2, default=[3, 5],
                        help="每属选取物种数范围 (min max)")
    parser.add_argument("--n-genera", type=int, default=0,
                        help="选取属的数量 (0=全部)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    by_genus = load_ref_info(args.ref_info)
    selected, log = select_per_genus(by_genus, args.n_per_genus[0], args.n_per_genus[1], args.n_genera, rng)
    write_output(selected, log, args.ref_fasta, args.outdir)


if __name__ == "__main__":
    main()
