#!/usr/bin/env python3
"""
构建病毒分类评估数据集 v3
  — 已知病毒 (--virus-dir) + 新病毒 (--novel-virus) → 覆盖度截取 + 突变
  — 新病毒 taxonomy 自动从 NCBI 补齐 → 输出 test_metadata_full.tsv

用法:
  python prep_build_class_eval_seqs.py \
      --virus-dir step1_eval_viruses/ \
      --ref-info final.cluster.ref_info.tsv \
      --ref-fasta final.cluster.ref.fasta \
      --novel-virus novel_viruses.fasta \
      --coverage-levels 100 90 80 70 60 50 40 --n-per-cov 2 \
      --mutation-rates 0.00 0.05 0.10 0.15 --n-per-mut 1 \
      --outdir step4_classification_eval/ --seed 42 \
      --email your@email.com
"""

import argparse, os, sys, random, csv, time, xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict, Counter
from urllib.request import urlopen
from urllib.parse import quote

from Bio import SeqIO
import pandas as pd

NUCLEOTIDES = ["A", "C", "G", "T"]
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
RANK_TO_ICTV = {
    'realm': 'Realm', 'kingdom': 'Kingdom', 'phylum': 'Phylum',
    'class': 'Class', 'order': 'Order', 'family': 'Family',
    'genus': 'Genus', 'species': 'Species',
}


# ══════════════════════════════════════
# 1. 加载
# ══════════════════════════════════════

def load_ref_info(path):
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
    print(f"[ref_info] {len(info)} records")
    return info


def load_viruses_from_dir(virus_dir, ref_info):
    """加载已知病毒"""
    files = sorted(Path(virus_dir).glob("*.fasta")) + sorted(Path(virus_dir).glob("*.fna"))
    seqs = []
    for f in files:
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            info = ref_info.get(acc, {})
            sp, g, fm = info.get("species", ""), info.get("genus", ""), info.get("family", "")
            if sp and g and fm:
                seqs.append((acc, str(rec.seq), sp, g, fm, len(rec.seq), "known"))
    print(f"[known] {len(seqs)} viruses from {virus_dir}")
    return seqs


def load_novel_viruses(fasta_path, ref_info, email):
    """加载新病毒 FASTA，缺失 taxonomy 的从 NCBI 在线补全"""
    seqs_raw = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        acc = rec.id.split()[0]
        info = ref_info.get(acc, {})
        sp = info.get("species", "")
        g = info.get("genus", "")
        fm = info.get("family", "")
        seqs_raw.append((acc, str(rec.seq), sp, g, fm, len(rec.seq)))

    # 找出需要 NCBI 查询的
    need_ncbi = [(acc, seq, sp, g, fm, ln)
                 for acc, seq, sp, g, fm, ln in seqs_raw
                 if not sp or not g or not fm]
    has_tax = [(acc, seq, sp, g, fm, ln, "novel")
               for acc, seq, sp, g, fm, ln in seqs_raw
               if sp and g and fm]

    if need_ncbi:
        print(f"[novel] {len(need_ncbi)}/{len(seqs_raw)} viruses need NCBI taxonomy query")
        filled = fetch_ncbi_taxonomy_batch(need_ncbi, email)
        has_tax.extend(filled)

    print(f"[novel] {len(has_tax)} viruses loaded")
    return has_tax


def fetch_ncbi_taxonomy_batch(novel_seqs, email):
    """批量从 NCBI 获取 taxonomy（带 500ms 延迟）"""
    results = []
    for i, (acc, seq, sp, g, fm, ln) in enumerate(novel_seqs):
        try:
            # Step 1: nuccore → taxid
            url = f"{NCBI_BASE}/efetch.fcgi?db=nuccore&id={quote(acc)}&rettype=xml&retmode=xml&email={quote(email)}"
            resp = urlopen(url, timeout=30)
            root = ET.fromstring(resp.read())

            # 提取 taxid
            taxid_el = root.find('.//GBSeq_taxid')
            taxid = taxid_el.text if taxid_el is not None and taxid_el.text != '0' else None
            if not taxid:
                for el in root.iter('GBQualifier'):
                    n_el, v_el = el.find('GBQualifier_name'), el.find('GBQualifier_value')
                    if n_el is not None and v_el is not None:
                        if n_el.text == 'db_xref' and v_el.text and v_el.text.startswith('taxon:'):
                            taxid = v_el.text.split(':')[1]
                            break
            if not taxid:
                print(f"  [{i+1}/{len(novel_seqs)}] {acc} FAILED (no taxid)")
                results.append((acc, seq, sp or f"novel_sp_{acc}",
                               g or "novel_genus", fm or "novel_family", ln, "novel"))
                continue

            # 提取 nuccore 中的 species/genus/family
            if not sp:
                org_el = root.find('.//GBSeq_organism')
                if org_el is not None and org_el.text:
                    sp = org_el.text.strip()
            if not fm:
                tax_el = root.find('.//GBSeq_taxonomy')
                if tax_el is not None and tax_el.text:
                    parts = [p.strip() for p in tax_el.text.split(';')]
                    # 倒数第二、第三个通常是 genus/family
                    if len(parts) >= 2 and not g:
                        g = parts[-2]
                    if len(parts) >= 3 and not fm:
                        fm = parts[-3]

            # Step 2: taxonomy → ranked lineage（高等级）
            time.sleep(0.4)
            url2 = f"{NCBI_BASE}/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml&email={quote(email)}"
            resp2 = urlopen(url2, timeout=30)
            root2 = ET.fromstring(resp2.read())
            taxon = root2.find('.//Taxon')

            high_tax = {'Realm': '', 'Kingdom': '', 'Phylum': '', 'Class': '', 'Order': ''}
            if taxon is not None:
                lex = taxon.find('LineageEx')
                if lex is not None:
                    for child in lex:
                        rk, nm = child.find('Rank'), child.find('ScientificName')
                        if rk is not None and nm is not None:
                            rank = rk.text.lower()
                            name = nm.text
                            if rank in RANK_TO_ICTV and RANK_TO_ICTV[rank] in high_tax:
                                high_tax[RANK_TO_ICTV[rank]] = name
                rk, nm = taxon.find('Rank'), taxon.find('ScientificName')
                if rk is not None and nm is not None:
                    rank = rk.text.lower()
                    if rank in RANK_TO_ICTV and rank == 'family' and not fm:
                        fm = nm.text
                    if rank in RANK_TO_ICTV and rank == 'genus' and not g:
                        g = nm.text
                    if rank in RANK_TO_ICTV and RANK_TO_ICTV[rank] in high_tax:
                        high_tax[RANK_TO_ICTV[rank]] = nm.text

            sp = sp or f"novel_sp_{acc}"
            g = g or "novel_genus"
            fm = fm or "novel_family"
            results.append((acc, seq, sp, g, fm, ln, "novel", high_tax))
            print(f"  [{i+1}/{len(novel_seqs)}] {acc} OK "
                  f"(sp={sp[:40]}, fm={fm}, Realm={high_tax.get('Realm', '')})")

        except Exception as e:
            print(f"  [{i+1}/{len(novel_seqs)}] {acc} FAILED ({e})")
            results.append((acc, seq, sp or f"novel_sp_{acc}",
                           g or "novel_genus", fm or "novel_family", ln, "novel"))

    return results


# ══════════════════════════════════════
# 2. 测试序列生成
# ══════════════════════════════════════

def mutate_sequence(seq_str, rate, rng):
    if rate == 0.0:
        return seq_str
    bases = list(seq_str)
    for i in range(len(bases)):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in NUCLEOTIDES if b != bases[i]])
    return "".join(bases)


def generate_all_tests(virus_seqs, coverage_levels, n_per_cov,
                       mutation_rates, n_per_mut, rng):
    """先突变，再从突变序列截取。生成：coverage + mut{rate} + mut{rate}_cov{cov}"""
    test_records, test_seqs = [], []
    idx = 0

    N_COV = len(coverage_levels) * n_per_cov
    N_MUT = len(mutation_rates) * n_per_mut
    N_MUT_COV = len(mutation_rates) * n_per_mut * len(coverage_levels) * n_per_cov

    for item in virus_seqs:
        acc, seq, sp, g, fm, ln, category = item[:7]
        extra_tax = item[7] if len(item) > 7 and isinstance(item[7], dict) else {}

        # ── 1. 原序列覆盖度截取 (type=coverage) ──
        for cov_pct in coverage_levels:
            for _ in range(n_per_cov):
                frag_len = max(500, int(ln * cov_pct / 100))
                frag_len = min(frag_len, ln)
                max_start = ln - frag_len
                start = rng.randint(0, max(1, max_start)) if max_start > 0 else 0
                frag = seq[start:start + frag_len]
                seq_id = f"test|cov{cov_pct}|{idx:04d}|src={acc}"
                test_seqs.append((seq_id, frag))
                rec = {
                    "seq_id": seq_id, "source_accession": acc,
                    "species": sp, "genus": g, "family": fm,
                    "coverage_pct": cov_pct, "mutation_rate_pct": 0,
                    "full_length": ln, "frag_length": len(frag),
                    "virus_type": category,
                    "mut_type": "mut0",
                    "cov_type": f"cov{cov_pct}",
                }
                rec.update(extra_tax)
                test_records.append(rec)
                idx += 1

        # ── 2. 突变 + 突变后覆盖度截取 ──
        for rate in mutation_rates:
            rate_pct = int(rate * 100)
            mut_type = f"mut{rate_pct}"

            for _ in range(n_per_mut):
                mutated = mutate_sequence(seq, rate, rng)

                # 2a. 突变全长
                seq_id = f"test|{mut_type}|{idx:04d}|src={acc}"
                test_seqs.append((seq_id, mutated))
                rec = {
                    "seq_id": seq_id, "source_accession": acc,
                    "species": sp, "genus": g, "family": fm,
                    "coverage_pct": 100, "mutation_rate_pct": rate_pct,
                    "full_length": ln, "frag_length": len(mutated),
                    "virus_type": category,
                    "mut_type": mut_type, "cov_type": "cov100",
                }
                rec.update(extra_tax)
                test_records.append(rec)
                idx += 1

                # 2b. 突变序列覆盖度截取
                for cov_pct in coverage_levels:
                    cov_type = f"cov{cov_pct}"
                    for _ in range(n_per_cov):
                        frag_len = max(500, int(ln * cov_pct / 100))
                        frag_len = min(frag_len, ln)
                        max_start = ln - frag_len
                        start = rng.randint(0, max(1, max_start)) if max_start > 0 else 0
                        frag = mutated[start:start + frag_len]
                        seq_id = f"test|{mut_type}_{cov_type}|{idx:04d}|src={acc}"
                        test_seqs.append((seq_id, frag))
                        rec = {
                            "seq_id": seq_id, "source_accession": acc,
                            "species": sp, "genus": g, "family": fm,
                            "coverage_pct": cov_pct,
                            "mutation_rate_pct": rate_pct,
                            "full_length": ln, "frag_length": len(frag),
                            "virus_type": category,
                            "mut_type": mut_type, "cov_type": cov_type,
                        }
                        rec.update(extra_tax)
                        test_records.append(rec)
                        idx += 1

    n_v = len(virus_seqs)
    total_expected = n_v * (N_COV + N_MUT + N_MUT_COV)
    print(f"[test] {len(test_records)} sequences ({n_v} viruses × "
          f"[{N_COV} cov + {N_MUT} mut + {N_MUT_COV} mut×cov])")
    return test_records, test_seqs


# ══════════════════════════════════════
# 3. 本地补全已知病毒高等级分类
# ══════════════════════════════════════

def enrich_known_taxonomy_local(meta_df, pred_tsv=None):
    """对 virus_type=known 的序列，用预测文件补全 Realm→Order"""
    if pred_tsv and os.path.exists(pred_tsv):
        pred = pd.read_csv(pred_tsv, sep='\t')
        sub = pred[pred['Family'].notna() & (pred['Family'] != '')]
        family_map = {}
        for family, grp in sub.groupby('Family'):
            for lvl in ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order']:
                vals = grp[lvl].dropna()
                vals = vals[vals != '']
                if len(vals) > 0:
                    family_map.setdefault(family, {})[lvl] = vals.mode().iloc[0]

        known_mask = meta_df['virus_type'] == 'known'
        for lvl in ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order']:
            meta_df.loc[known_mask, lvl] = meta_df.loc[known_mask, 'family'].map(
                lambda f: family_map.get(f, {}).get(lvl, ''))
        print(f"[enrich-known] {len(family_map)} families mapped for known viruses")
    else:
        for lvl in ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order']:
            if lvl not in meta_df.columns:
                meta_df[lvl] = ''
    return meta_df


def enrich_taxonomy_ncbi(meta_df, email, cache_path=None):
    """NCBI Entrez 两阶段查询已知病毒的高等级分类"""
    cache = {}
    if cache_path and os.path.exists(cache_path):
        existing = pd.read_csv(cache_path, sep='\t')
        for _, row in existing.iterrows():
            cache[row['accession']] = row.to_dict()

    known_mask = meta_df['virus_type'] == 'known'
    accessions = meta_df.loc[known_mask, 'source_accession'].dropna().unique()
    new_acc = [a for a in accessions if a not in cache]

    if not new_acc:
        print("[enrich-ncbi-known] All cached, skipping")
    else:
        print(f"[enrich-ncbi-known] Querying {len(new_acc)} accessions...")
        for i, acc in enumerate(new_acc):
            try:
                url = (f"{NCBI_BASE}/efetch.fcgi?db=nuccore&id={quote(acc)}"
                       f"&rettype=xml&retmode=xml&email={quote(email)}")
                resp = urlopen(url, timeout=30)
                root = ET.fromstring(resp.read())
                taxid = None
                for el in root.iter('GBSeq_taxid'):
                    if el.text and el.text != '0':
                        taxid = el.text
                        break
                if not taxid:
                    for el in root.iter('GBQualifier'):
                        n_el, v_el = el.find('GBQualifier_name'), el.find('GBQualifier_value')
                        if n_el is not None and v_el is not None:
                            if n_el.text == 'db_xref' and v_el.text and v_el.text.startswith('taxon:'):
                                taxid = v_el.text.split(':')[1]
                                break
                if not taxid:
                    raise ValueError("no taxid")

                time.sleep(0.35)
                url2 = f"{NCBI_BASE}/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml&email={quote(email)}"
                resp2 = urlopen(url2, timeout=30)
                root2 = ET.fromstring(resp2.read())
                taxon = root2.find('.//Taxon')
                entry = {'accession': acc, 'Realm': '', 'Kingdom': '',
                         'Phylum': '', 'Class': '', 'Order': ''}
                if taxon is not None:
                    lex = taxon.find('LineageEx')
                    if lex is not None:
                        for child in lex:
                            rk, nm = child.find('Rank'), child.find('ScientificName')
                            if rk is not None and nm is not None:
                                rank = rk.text.lower()
                                if rank in RANK_TO_ICTV and RANK_TO_ICTV[rank] in entry:
                                    entry[RANK_TO_ICTV[rank]] = nm.text
                    rk, nm = taxon.find('Rank'), taxon.find('ScientificName')
                    if rk is not None and nm is not None:
                        rank = rk.text.lower()
                        if rank in RANK_TO_ICTV and RANK_TO_ICTV[rank] in entry:
                            entry[RANK_TO_ICTV[rank]] = nm.text
                cache[acc] = entry
                print(f"  [{i+1}/{len(new_acc)}] {acc} OK (Realm={entry['Realm']})")
            except Exception as e:
                print(f"  [{i+1}/{len(new_acc)}] {acc} FAILED ({e})")
                cache[acc] = {'accession': acc, 'Realm': '', 'Kingdom': '',
                              'Phylum': '', 'Class': '', 'Order': ''}

        if cache_path:
            pd.DataFrame(list(cache.values())).to_csv(cache_path, sep='\t', index=False)
            print(f"[enrich-ncbi-known] Cache: {cache_path}")

    for lvl in ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order']:
        meta_df.loc[known_mask, lvl] = meta_df.loc[known_mask, 'source_accession'].map(
            lambda a: cache.get(a, {}).get(lvl, ''))
    return meta_df


def ensure_tax_columns(meta_df):
    """确保所有 taxonomy 列存在"""
    for lvl in ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order']:
        if lvl not in meta_df.columns:
            meta_df[lvl] = ''
    return meta_df


# ══════════════════════════════════════
# 4. 假阳性（decoy）序列
# ══════════════════════════════════════

def load_decoys(conserved_fasta, eve_fasta, host_fasta, n_decoys, frag_lengths, rng):
    """加载假阳性序列：>500bp，随机截取片段匹配阳性长度分布，不突变"""
    decoy_records, decoy_seqs = [], []
    idx = 0
    # 参考 prep_master_eval_dataset.py 的 D/E/H 阴性组
    sources = [
        ("pfam",  conserved_fasta, True),   # 保守结构域目录（多文件）
        ("eve",   eve_fasta,      False),
        ("host",  host_fasta,     False),
    ]

    for virus_type, fasta_path, is_dir in sources:
        if not fasta_path or not os.path.exists(fasta_path):
            print(f"[decoy] {virus_type}: file not found, skipping")
            continue

        # 加载 >500bp 的序列
        raw_seqs = []
        if is_dir and os.path.isdir(fasta_path):
            for f in sorted(Path(fasta_path).glob("*.fasta")):
                for rec in SeqIO.parse(f, "fasta"):
                    if len(rec.seq) > 500:
                        raw_seqs.append((rec.id.split()[0], str(rec.seq)))
        elif os.path.isfile(fasta_path):
            for rec in SeqIO.parse(fasta_path, "fasta"):
                if len(rec.seq) > 500:
                    raw_seqs.append((rec.id.split()[0], str(rec.seq)))
        else:
            print(f"[decoy] {virus_type}: '{fasta_path}' is not a valid file/dir")
            continue

        if len(raw_seqs) == 0:
            print(f"[decoy] {virus_type}: no sequences >500bp, skipping")
            continue

        n = min(n_decoys, len(raw_seqs))
        sampled = rng.sample(raw_seqs, n) if n < len(raw_seqs) else raw_seqs

        for acc, full_seq in sampled:
            # 随机截取片段，长度匹配阳性序列分布（参考 master 脚本）
            flen = min(rng.choice(frag_lengths), len(full_seq))
            flen = max(500, flen)
            start = rng.randint(0, max(1, len(full_seq) - flen))
            frag = full_seq[start:start + flen]

            seq_id = f"test|decoy_{virus_type}|{idx:04d}|src={acc}"
            decoy_seqs.append((seq_id, frag))
            decoy_records.append({
                "seq_id": seq_id,
                "source_accession": acc,
                "species": f"decoy_{virus_type}",
                "genus": f"decoy_{virus_type}",
                "family": f"decoy_{virus_type}",
                "coverage_pct": 100,
                "mutation_rate_pct": 0,
                "full_length": len(full_seq),
                "frag_length": len(frag),
                "virus_type": virus_type,
                "mut_type": "",
                "cov_type": "",
                "Realm": "", "Kingdom": "", "Phylum": "", "Class": "", "Order": "",
            })
            idx += 1
        print(f"[decoy] {virus_type}: {n} sequences")

    print(f"[decoy] Total: {len(decoy_records)} decoy sequences")
    return decoy_records, decoy_seqs


# ══════════════════════════════════════
# 5. 主函数
# ══════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="构建病毒分类评估数据集 v3")
    parser.add_argument("--virus-dir", required=True, help="已知病毒 FASTA 目录")
    parser.add_argument("--novel-virus", default=None,
                        help="新病毒 FASTA 文件（单个文件）")
    parser.add_argument("--ref-info", required=True)
    parser.add_argument("--ref-fasta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument("--coverage-levels", type=int, nargs="+",
                        default=[100, 90, 80, 70, 60, 50, 40])
    parser.add_argument("--n-per-cov", type=int, default=2)

    parser.add_argument("--mutation-rates", type=float, nargs="+",
                        default=[0.00, 0.05, 0.10, 0.15])
    parser.add_argument("--n-per-mut", type=int, default=1)

    parser.add_argument("--email", default="user@example.com",
                        help="NCBI email（查询新病毒 taxonomy 时需要）")

    parser.add_argument("--conserved-fasta", default=None,
                        help="保守结构域/Pfam FASTA 目录 → virus_type=pfam")
    parser.add_argument("--eve-fasta", default=None,
                        help="EVE 内源性病毒元件 FASTA → virus_type=eve")
    parser.add_argument("--host-fasta", default=None,
                        help="宿主基因组 FASTA → virus_type=host")
    parser.add_argument("--n-decoys", type=int, default=300,
                        help="每类假阳性序列数量")

    # 已知病毒的高等级补全（本地模式，需要预测文件）
    parser.add_argument("--enrich-method", choices=["none", "local", "ncbi"],
                        default="ncbi", help="已知病毒的高等级分类补全方式")
    parser.add_argument("--enrich-pred", default=None,
                        help="整合分类 TSV（--enrich-method local 时使用）")
    parser.add_argument("--enrich-cache", default=None,
                        help="NCBI 缓存文件（--enrich-method ncbi 时使用）")

    args = parser.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # 1. 加载
    ref_info = load_ref_info(args.ref_info)

    all_virus_seqs = load_viruses_from_dir(args.virus_dir, ref_info)

    if args.novel_virus and os.path.exists(args.novel_virus):
        novel_seqs = load_novel_viruses(args.novel_virus, ref_info, args.email)
        all_virus_seqs.extend(novel_seqs)

    if len(all_virus_seqs) < 3:
        print("[ERROR] Too few virus sequences!")
        sys.exit(1)

    n_known = sum(1 for s in all_virus_seqs if s[6] == "known")
    n_novel = sum(1 for s in all_virus_seqs if s[6] == "novel")
    print(f"[total] {len(all_virus_seqs)} viruses (known={n_known}, novel={n_novel})")

    # 2. 生成测试序列（已知+新病毒：覆盖度截取+突变）
    test_records, test_seqs = generate_all_tests(
        all_virus_seqs, args.coverage_levels, args.n_per_cov,
        args.mutation_rates, args.n_per_mut, rng)

    # 2b. 加载假阳性序列（>500bp，随机截取匹配阳性长度分布）
    frag_lengths = [r["frag_length"] for r in test_records] if test_records else [500]
    decoy_records, decoy_seqs = load_decoys(
        args.conserved_fasta, args.eve_fasta, args.host_fasta,
        args.n_decoys, frag_lengths, rng)
    test_records.extend(decoy_records)
    test_seqs.extend(decoy_seqs)

    # 3. 写 FASTAs
    test_dir = os.path.join(args.outdir, "test_sequences")
    os.makedirs(test_dir, exist_ok=True)
    for seq_id, seq_str in test_seqs:
        fname = seq_id.replace("|", "_") + ".fasta"
        with open(os.path.join(test_dir, fname), "w") as f:
            f.write(f">{seq_id}\n{seq_str}\n")

    merged = os.path.join(args.outdir, "test_sequences_merged.fasta")
    with open(merged, "w") as f:
        for seq_id, seq_str in test_seqs:
            f.write(f">{seq_id}\n{seq_str}\n")

    # 4. 元数据 + 补全分类
    meta_df = pd.DataFrame(test_records)
    meta_df = ensure_tax_columns(meta_df)
    if args.enrich_method == "local":
        meta_df = enrich_known_taxonomy_local(meta_df, args.enrich_pred)
    elif args.enrich_method == "ncbi":
        meta_df = enrich_taxonomy_ncbi(meta_df, args.email, args.enrich_cache)
    # else "none": 保持空列

    meta_path = os.path.join(args.outdir, "test_metadata_full.tsv")
    meta_df.to_csv(meta_path, sep="\t", index=False)
    print(f"[meta] {meta_path} ({len(meta_df)} rows, {len(meta_df.columns)} cols)")

    # 5. 去泄漏参考库
    exclude_species = set(r["species"] for r in test_records)
    if os.path.isdir(args.virus_dir):
        db_files = sorted(Path(args.virus_dir).glob("*.fasta")) + sorted(Path(args.virus_dir).glob("*.fna"))
    else:
        db_files = [Path(args.virus_dir)]
    kept = []
    for f in db_files:
        for rec in SeqIO.parse(f, "fasta"):
            acc = rec.id.split()[0]
            info = ref_info.get(acc, {})
            if info.get("species", "") in exclude_species:
                continue
            kept.append(rec)

    db_path = os.path.join(args.outdir, "db_sequences.fasta")
    SeqIO.write(kept, db_path, "fasta")
    print(f"[db] {len(kept)} sequences (excluded species={len(exclude_species)})")

    # 6. 统计
    print(f"\n[DONE] {len(test_records)} test sequences")
    for vt in ["known", "novel"]:
        recs_vt = [r for r in test_records if r["virus_type"] == vt]
        if not recs_vt:
            continue
        muts = sorted(set(r["mut_type"] for r in recs_vt))
        covs = sorted(set(r["cov_type"] for r in recs_vt),
                       key=lambda x: -int(x.replace("cov", "")))
        n_full = sum(1 for r in recs_vt if r.get("coverage_pct", 100) == 100)
        n_trunc = sum(1 for r in recs_vt if r.get("coverage_pct", 100) < 100)
        print(f"  {vt:6s}  {len(recs_vt):5d} seqs  "
              f"mut={muts}  cov={covs}  (full={n_full}, truncated={n_trunc})")
    for vt in ["pfam", "eve", "host"]:
        n = sum(1 for r in test_records if r["virus_type"] == vt)
        if n > 0:
            print(f"  {vt:6s}  {n:5d} seqs  (decoy, no mutation/coverage)")
    print(f"  Output: {args.outdir}")


if __name__ == "__main__":
    main()
