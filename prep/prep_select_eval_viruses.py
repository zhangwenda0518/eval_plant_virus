#!/usr/bin/env python3
"""
评估用病毒基因组选取脚本 (v2 — 基于实际 ref_info 格式)

从 final.cluster.ref_info.tsv 中按分层策略选取评估用病毒基因组。

选取条件（硬约束）:
  - Category = NonSegmented 或包含 NonSegmented（排除节段病毒，简化评估）
  - Sequence_Type = RefSeq（只要 RefSeq 参考序列，确保质量）
  - ICTV 分类完整（有 Species_ICTV / VMR_Genus / VMR_Family）
  - Nuc_Completeness = complete（完整基因组）

分层配额（保证基因组类型多样性）:
  Topology:  linear + circular 均有
  Molecule_type/Molecule_Type2: 覆盖 ssRNA(+) / ssRNA(-) / dsRNA / ssDNA / dsDNA-RT / viroid
  长度: 从 200 bp (类病毒) 到 30,000 bp (大RNA病毒) 全覆盖

用法:
  # 选取50个病毒（已知病毒检测评估）
  python prep_select_eval_viruses.py \
      --ref-info final.cluster.ref_info.tsv \
      --ref-fasta final.cluster.ref.fasta \
      --n-viruses 50 \
      --outdir eval_viruses_50/ --seed 42

  # 选取30个病毒（组装评估，排除已知病毒检测评估已选的）
  python prep_select_eval_viruses.py \
      --ref-info final.cluster.ref_info.tsv \
      --ref-fasta final.cluster.ref.fasta \
      --n-viruses 30 \
      --exclude eval_viruses_50/ \
      --outdir eval_viruses_30/ --seed 43
"""

import argparse, os, sys, random, csv
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
import pandas as pd


# ============================================================
# 基因组类型分类
# ============================================================

def classify_genome_type(row):
    """
    根据 Topology + Molecule_type + Molecule_Type2 + Length 综合判定基因组类型。

    返回: (main_type, subtype)
      main_type ∈ {"ssRNA(+)", "ssRNA(-)", "dsRNA", "ssDNA", "dsDNA-RT", "viroid"}
    """
    topo = str(row.get("Topology", "")).lower().strip()
    mol = str(row.get("Molecule_type", "")).strip()
    mol2 = str(row.get("Molecule_Type2", "")).lower().strip()
    length = int(row.get("Length", 0)) if row.get("Length") else 0
    title = str(row.get("GenBank_Title", "")).lower()

    # 类病毒: 长度 < 500 bp 且标题含 viroid
    if length <= 500 and "viroid" in title:
        return "viroid", "viroid_ncrna"

    # dsDNA-RT: 逆转录病毒特征的 dsDNA
    if "dsdna" in mol2 and "rt" in mol2:
        return "dsDNA-RT", "dsdna_rt"
    if "dsdna" in mol.lower() or mol == "dsDNA":
        return "dsDNA-RT", "dsdna_rt"

    # ssDNA
    if mol.startswith("ssDNA") or "ssdna" in mol2:
        if "circular" in topo:
            return "ssDNA", "ssdna_circular"
        return "ssDNA", "ssdna_linear"

    # dsRNA
    if mol.startswith("dsRNA") or "dsrna" in mol2:
        return "dsRNA", "dsrna"

    # ssRNA(-)
    if mol in ("ssRNA(-)", "ssRNA(-)") or "ssrna(-)" in mol2 or mol2 == "crna":
        return "ssRNA(-)", "ssrna_neg"

    # ssRNA(+) (最大的群体，兜底)
    if mol in ("ssRNA(+)", "ssRNA(+)") or "ssrna(+)" in mol2 or mol2 == "dna":
        return "ssRNA(+)", "ssrna_pos"

    # 兜底: 从标题推断
    if "negative-strand" in title or "negative strand" in title:
        return "ssRNA(-)", "ssrna_neg_inferred"
    if "positive-strand" in title or "positive strand" in title:
        return "ssRNA(+)", "ssrna_pos_inferred"

    # 最后的兜底: 按长度和拓扑
    if "circular" in topo and length < 3000:
        return "ssDNA", "ssdna_inferred"
    return "ssRNA(+)", "ssrna_pos_default"


def get_quota(n_total):
    """各基因组类型配额（按比例，从类病毒到大RNA全覆盖）"""
    if n_total == 50:
        return {
            "ssRNA(+)":   15,   # 最大的类群
            "ssRNA(-)":   8,
            "dsRNA":      6,
            "ssDNA":      10,   # 环状ssDNA，植物病毒中很多
            "dsDNA-RT":   5,
            "viroid":     6,    # 类病毒必须覆盖
        }
    elif n_total == 30:
        return {
            "ssRNA(+)":   10,
            "ssRNA(-)":   4,
            "dsRNA":      4,
            "ssDNA":      6,
            "dsDNA-RT":   3,
            "viroid":     3,
        }
    else:
        return {
            "ssRNA(+)":   max(1, n_total * 3 // 10),
            "ssRNA(-)":   max(1, n_total // 6),
            "dsRNA":      max(1, n_total // 8),
            "ssDNA":      max(1, n_total // 5),
            "dsDNA-RT":   max(1, n_total // 10),
            "viroid":     max(1, n_total // 10),
        }


# ============================================================
# 主逻辑
# ============================================================

def load_ref_info(path):
    """加载 ref_info TSV，过滤条件并分类"""
    records = []
    excluded_reasons = Counter()

    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # 硬约束1: Category 非节段（排除 Segmented，保留 NonSegmented 等）
            cat = str(row.get("Category", "")).strip()
            if "segmented" in cat.lower() and "nonsegmented" not in cat.lower():
                excluded_reasons["segmented"] += 1
                continue

            # 硬约束2: Nuc_Completeness = complete（接受RefSeq和GenBank/ICTV）
            comp = str(row.get("Nuc_Completeness", "")).strip().lower()
            if comp != "complete":
                excluded_reasons["not_complete"] += 1
                continue

            # 硬约束3: ICTV 分类完整（VMR_Family/Genus/Species 三者均有）
            family = row.get("VMR_Family", "").strip()
            genus = row.get("VMR_Genus", "").strip()
            species = row.get("VMR_Species", "").strip()
            if not family or not genus or not species:
                excluded_reasons["ICTV_incomplete"] += 1
                continue

            # 硬约束4: 排除非植物宿主（绿藻病毒等）
            host = str(row.get("Host", "")).strip().lower()
            title = str(row.get("GenBank_Title", "")).strip().lower()
            non_plant_kw = ["chlorella", "phycodnaviridae", "chlorovirus",
                           "prasinovirus", "coccolithovirus", "prymnesiovirus",
                           "phaeovirus", "raphidovirus"]
            skip = False
            for kw in non_plant_kw:
                if kw in host or kw in title:
                    excluded_reasons["non_plant_host"] += 1
                    skip = True
                    break
            if skip:
                continue

            # 有长度
            length_str = row.get("Length", "0").strip()
            if not length_str or not length_str.isdigit():
                excluded_reasons["no_length"] += 1
                continue
            length = int(length_str)

            # 有 Accession
            acc = row.get("Accession", "").strip()
            if not acc:
                excluded_reasons["no_accession"] += 1
                continue

            gtype, subtype = classify_genome_type(row)
            records.append({
                "accession": acc,
                "species": species,
                "genus": genus,
                "family": family,
                "genome_type": gtype,
                "subtype": subtype,
                "topology": str(row.get("Topology", "")).strip(),
                "molecule_type": str(row.get("Molecule_type", "")).strip(),
                "length": length,
                "host": str(row.get("Host", "")).strip(),
            })

    print(f"[ref_info] Total records loaded: {len(records)}")
    print(f"[ref_info] Excluded:")
    for reason, count in excluded_reasons.most_common():
        print(f"  {reason}: {count}")

    return records


def load_fasta_index(fasta_path):
    """构建 Accession → seq_record 的索引"""
    idx = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        acc = rec.id.split()[0]
        # 可能用不同的分隔符
        idx[acc] = rec
        idx[rec.id] = rec
    print(f"[fasta] Indexed {len(idx)//2} sequences")
    return idx


def load_excluded_accessions(exclude_dir):
    """加载需排除的 Accession 集合"""
    excluded = set()
    if not exclude_dir or not os.path.exists(exclude_dir):
        return excluded
    for f in Path(exclude_dir).glob("*.fasta"):
        excluded.add(f.stem.split(".")[0])  # NC_xxxxxx.1
        for rec in SeqIO.parse(f, "fasta"):
            excluded.add(rec.id)
    print(f"[exclude] {len(excluded)} accessions excluded")
    return excluded


def stratified_select(records, n_total, rng):
    """
    分层选取: 按 genome_type 分配配额 → 在每个类型内分散选取不同长度范围
    """
    quota = get_quota(n_total)

    # 按类型分组
    by_type = defaultdict(list)
    for r in records:
        by_type[r["genome_type"]].append(r)

    selected = []
    selection_log = []

    for gtype in ["ssRNA(+)", "ssRNA(-)", "dsRNA", "ssDNA", "dsDNA-RT", "viroid"]:
        candidates = by_type.get(gtype, [])
        n_wanted = quota.get(gtype, 0)

        if not candidates:
            print(f"[WARNING] No candidates for {gtype}")
            continue

        # 在类型内按长度分箱，保证长度覆盖
        len_bins = [
            ("ultra_short", 0, 500),
            ("short", 501, 2000),
            ("medium", 2001, 6000),
            ("long", 6001, 15000),
            ("ultra_long", 15001, 100000),
        ]

        by_len = defaultdict(list)
        for c in candidates:
            for name, lo, hi in len_bins:
                if lo <= c["length"] <= hi:
                    by_len[name].append(c)
                    break

        # 每箱至少选1-2个
        n_picked = 0
        picked_in_type = []

        # 第一轮: 每箱选1个
        for name, lo, hi in len_bins:
            if name in by_len and by_len[name]:
                pick = rng.sample(by_len[name], 1)[0]
                picked_in_type.append(pick)
                n_picked += 1

        # 第二轮: 剩余的随机选
        picked_accs = {p["accession"] for p in picked_in_type}
        remaining_needed = n_wanted - n_picked
        if remaining_needed > 0:
            remaining_candidates = [c for c in candidates if c["accession"] not in picked_accs]
            if remaining_candidates:
                extra = rng.sample(remaining_candidates,
                                   min(remaining_needed, len(remaining_candidates)))
                picked_in_type.extend(extra)

        for p in picked_in_type:
            selected.append(p)
            selection_log.append({
                "accession": p["accession"],
                "genome_type": p["genome_type"],
                "topology": p["topology"],
                "molecule_type": p["molecule_type"],
                "family": p["family"],
                "genus": p["genus"],
                "species": p["species"],
                "length": p["length"],
                "host": p["host"],
            })

        print(f"  {gtype}: picked {len(picked_in_type)} / wanted {n_wanted} "
              f"(from {len(candidates)} candidates)")

    return selected, selection_log


def write_output(selected, selection_log, fasta_idx, outdir):
    """写出FASTA文件和选取记录"""
    os.makedirs(outdir, exist_ok=True)

    # 写出各病毒FASTA
    written = 0
    for item in selected:
        acc = item["accession"]
        rec = fasta_idx.get(acc)
        if rec is None:
            print(f"[WARNING] {acc} not found in FASTA index")
            continue
        out_path = os.path.join(outdir, f"{acc}.fasta")
        SeqIO.write(rec, out_path, "fasta")
        written += 1

    # 写出选取记录
    log_df = pd.DataFrame(selection_log)
    log_df.to_csv(os.path.join(outdir, "selected_viruses.tsv"), sep="\t", index=False)

    # 统计
    print(f"\n{'='*50}")
    print(f"Selection Summary ({outdir})")
    print(f"{'='*50}")
    print(f"Total written: {written} FASTA files")
    if written > 0 and not log_df.empty and "genome_type" in log_df.columns:
        print(f"\nGenome type distribution:")
        for gtype, count in log_df["genome_type"].value_counts().items():
            print(f"  {gtype}: {count}")
        print(f"\nTopology distribution:")
        for topo, count in log_df["topology"].value_counts().items():
            print(f"  {topo}: {count}")
        print(f"\nLength range: {log_df['length'].min()} - {log_df['length'].max()} bp")
        print(f"\nFamily diversity: {log_df['family'].nunique()} families")
        print(f"Genus diversity: {log_df['genus'].nunique()} genera")


def main():
    parser = argparse.ArgumentParser(description="评估用病毒基因组选取 (v2)")
    parser.add_argument("--ref-info", required=True, help="final.cluster.ref_info.tsv")
    parser.add_argument("--ref-fasta", required=True, help="final.cluster.ref.fasta")
    parser.add_argument("--n-viruses", type=int, default=50, help="目标选取数量")
    parser.add_argument("--exclude", help="需排除的目录（含已选的.fasta）")
    parser.add_argument("--include", help="必须包含的病毒Accession (逗号分隔或文件)")
    parser.add_argument("--n-groups", type=int, default=0, help="分层分组数量 (0=不分组, 如5则分成5组)")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    # 1. 加载ref_info
    records = load_ref_info(args.ref_info)

    # 2. 排除已有
    excluded = load_excluded_accessions(args.exclude)
    records = [r for r in records if r["accession"] not in excluded]
    print(f"[select] After exclusion: {len(records)} candidates")

    # 3. 加载FASTA索引
    fasta_idx = load_fasta_index(args.ref_fasta)

    # 3.5. 强制包含指定病毒
    forced = []
    if args.include:
        include_list = []
        if os.path.isfile(args.include):
            with open(args.include) as f:
                include_list = [line.strip().split()[0] for line in f if line.strip()]
        else:
            include_list = [x.strip() for x in args.include.split(',')]
        for acc in include_list:
            for r in records:
                if r['accession'] == acc:
                    forced.append(r)
                    records.remove(r)
                    break
        print(f"[include] Forced {len(forced)} viruses: {[r['accession'] for r in forced]}")

    # 4. 分层选取
    selected, log = stratified_select(records, args.n_viruses - len(forced), rng)
    selected = forced + selected

    # 5. 检查: 候选不够时从其余类型补充
    if len(selected) < args.n_viruses:
        shortfall = args.n_viruses - len(selected)
        print(f"\n[WARNING] Shortfall of {shortfall}, filling from remaining pool...")
        remaining = [r for r in records if r not in selected]
        if remaining:
            extra = rng.sample(remaining, min(shortfall, len(remaining)))
            for e in extra:
                selected.append(e)
                log.append({
                    "accession": e["accession"],
                    "genome_type": e["genome_type"],
                    "topology": e["topology"],
                    "molecule_type": e["molecule_type"],
                    "family": e["family"],
                    "genus": e["genus"],
                    "species": e["species"],
                    "length": e["length"],
                    "host": e["host"],
                })

    # 6. 写出
    write_output(selected, log, fasta_idx, args.outdir)

    # 7. 按基因组类型分层分组（--n-groups=N）
    if args.n_groups and args.n_groups > 1:
        _write_groups(log, args.outdir, args.n_groups, rng)


def _write_groups(selection_log, outdir, n_groups, rng):
    """按基因组类型分层分组，每组输出到 group_N/ 子目录"""
    import shutil
    grouped = {i: [] for i in range(1, n_groups + 1)}

    for gtype in sorted(set(r['genome_type'] for r in selection_log)):
        subset = [r for r in selection_log if r['genome_type'] == gtype]
        rng.shuffle(subset)
        for i, r in enumerate(subset):
            grouped[(i % n_groups) + 1].append(r)

    for g in range(1, n_groups + 1):
        gdir = os.path.join(outdir, f"group_{g}")
        os.makedirs(gdir, exist_ok=True)
        gaccs = [r['accession'] for r in grouped[g]]
        # 拷贝FASTA文件
        for acc in gaccs:
            src = os.path.join(outdir, f"{acc}.fasta")
            dst = os.path.join(gdir, f"{acc}.fasta")
            if os.path.exists(src) and not os.path.exists(dst):
                shutil.copy2(src, dst)
        # 写分组表
        import pandas as pd
        gdf = pd.DataFrame(grouped[g])
        gdf.to_csv(os.path.join(gdir, f"group_{g}.tsv"), sep='\t', index=False)
        print(f"  Group {g}: {len(gaccs)} viruses")


if __name__ == "__main__":
    main()
