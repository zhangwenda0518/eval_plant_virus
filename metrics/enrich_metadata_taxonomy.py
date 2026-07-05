#!/usr/bin/env python3
"""
用 NCBI Entrez 根据 source_accession 补全 metadata 的高等级分类
通过 taxonomy XML 获取带 rank 标记的完整 lineage，精确映射 ICTV 等级。

用法:
  python enrich_metadata_taxonomy.py \
      --meta test_metadata.tsv \
      --out test_metadata_full.tsv \
      --email your@email.com

本地方案:
  --local --pred final_integrated_classification.tsv
"""

import argparse, os, sys, time, xml.etree.ElementTree as ET
import pandas as pd
from urllib.request import urlopen
from urllib.parse import quote

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TARGET_LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order']

# NCBI rank → ICTV level (只取标准等级，跳过 sub/super/infra/no rank)
RANK_TO_ICTV = {
    'realm':       'Realm',
    'kingdom':     'Kingdom',
    'phylum':      'Phylum',
    'class':       'Class',
    'order':       'Order',
    'family':      'Family',
    'genus':       'Genus',
    'species':     'Species',
}


def get_taxid(accession, email):
    """从 nuccore XML 中提取 taxid"""
    url = (f"{NCBI_BASE}/efetch.fcgi?db=nuccore&id={quote(accession)}"
           f"&rettype=xml&retmode=xml&email={quote(email)}")
    try:
        resp = urlopen(url, timeout=30)
        root = ET.fromstring(resp.read())
        for el in root.iter('GBSeq_taxid'):
            if el.text and el.text != '0':
                return el.text
        # 兜底：从 GBSeq_feature-table 的 qualifier 里找 db_xref taxon
        for el in root.iter('GBQualifier'):
            name_el = el.find('GBQualifier_name')
            val_el = el.find('GBQualifier_value')
            if name_el is not None and val_el is not None:
                if name_el.text == 'db_xref' and val_el.text and val_el.text.startswith('taxon:'):
                    return val_el.text.split(':')[1]
    except Exception as e:
        print(f"  [WARN] get_taxid failed for {accession}: {e}")
    return None


def get_ranked_lineage(taxid, email):
    """从 taxonomy XML 获取带 rank 标记的 lineage: [(rank, name), ...]"""
    url = (f"{NCBI_BASE}/efetch.fcgi?db=taxonomy&id={taxid}"
           f"&retmode=xml&email={quote(email)}")
    try:
        resp = urlopen(url, timeout=30)
        root = ET.fromstring(resp.read())
        taxon = root.find('.//Taxon')
        if taxon is None:
            return None

        lineage = []
        # LineageEx 中的祖先节点
        lineage_ex = taxon.find('LineageEx')
        if lineage_ex is not None:
            for child in lineage_ex:
                rank_el = child.find('Rank')
                name_el = child.find('ScientificName')
                if rank_el is not None and name_el is not None:
                    lineage.append((rank_el.text.lower(), name_el.text))

        # 加上自身
        rank_el = taxon.find('Rank')
        name_el = taxon.find('ScientificName')
        if rank_el is not None and name_el is not None:
            lineage.append((rank_el.text.lower(), name_el.text))

        return lineage
    except Exception as e:
        print(f"  [WARN] get_ranked_lineage failed for taxid {taxid}: {e}")
    return None


def lineage_to_levels(ranked_lineage):
    """将带 rank 的 lineage 映射为 {Realm, Kingdom, Phylum, Class, Order}"""
    result = {lvl: '' for lvl in TARGET_LEVELS}
    if not ranked_lineage:
        return result
    for rank, name in ranked_lineage:
        if rank in RANK_TO_ICTV:
            lvl = RANK_TO_ICTV[rank]
            if lvl in result:
                result[lvl] = name
    return result


def main():
    parser = argparse.ArgumentParser(description='补全 metadata 高等级分类')
    parser.add_argument('--meta', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--email', default='user@example.com')
    parser.add_argument('--cache', default=None)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--pred', default=None)
    args = parser.parse_args()

    meta = pd.read_csv(args.meta, sep='\t')

    if args.local and args.pred:
        pred = pd.read_csv(args.pred, sep='\t')
        sub = pred[pred['Family'].notna() & (pred['Family'] != '')]
        family_map = {}
        for family, grp in sub.groupby('Family'):
            for lvl in TARGET_LEVELS:
                vals = grp[lvl].dropna()
                vals = vals[vals != '']
                if len(vals) == 0:
                    continue
                mode_val = vals.mode()
                if len(mode_val) > 0:
                    family_map.setdefault(family, {})[lvl] = mode_val.iloc[0]
        print(f"Built taxonomy map for {len(family_map)} families")
        meta_fams = set(meta['family'].dropna().unique())
        covered = meta_fams & set(family_map.keys())
        missing = meta_fams - set(family_map.keys())
        print(f"Covered: {len(covered)}, Missing: {sorted(missing) if missing else 'None'}")
        for lvl in TARGET_LEVELS:
            meta[lvl] = meta['family'].map(lambda f: family_map.get(f, {}).get(lvl, ''))
    else:
        cache = {}
        if args.cache and os.path.exists(args.cache):
            existing = pd.read_csv(args.cache, sep='\t')
            for _, row in existing.iterrows():
                cache[row['accession']] = row.to_dict()

        accessions = meta['source_accession'].dropna().unique()
        new_queries = [a for a in accessions if a not in cache]
        print(f"Total unique accessions: {len(accessions)}")
        print(f"Cached: {len(cache)}, Need to query: {len(new_queries)}")

        for i, acc in enumerate(new_queries):
            print(f"  [{i+1}/{len(new_queries)}] {acc} ...", end=' ', flush=True)

            # Step 1: nuccore → taxid
            taxid = get_taxid(acc, args.email)
            if not taxid:
                print("FAILED (no taxid)")
                cache[acc] = {'accession': acc, 'Realm': '', 'Kingdom': '',
                              'Phylum': '', 'Class': '', 'Order': ''}
                continue

            # Step 2: taxonomy → ranked lineage
            time.sleep(0.35)
            ranked = get_ranked_lineage(taxid, args.email)
            if ranked:
                cache[acc] = lineage_to_levels(ranked)
                cache[acc]['accession'] = acc
                realm = cache[acc].get('Realm', '')
                order = cache[acc].get('Order', '')
                print(f"OK (Realm={realm}, Order={order})")
            else:
                print("FAILED (no lineage)")
                cache[acc] = {'accession': acc, 'Realm': '', 'Kingdom': '',
                              'Phylum': '', 'Class': '', 'Order': ''}

        if args.cache:
            cache_df = pd.DataFrame(list(cache.values()))
            cache_df.to_csv(args.cache, sep='\t', index=False)
            print(f"Cache saved: {args.cache}")

        lookup = {a: {lvl: cache[a].get(lvl, '') for lvl in TARGET_LEVELS} for a in cache}
        for lvl in TARGET_LEVELS:
            meta[lvl] = meta['source_accession'].map(
                lambda a: lookup.get(a, {}).get(lvl, ''))

    meta.to_csv(args.out, sep='\t', index=False)
    print(f"\nDone: {args.out}")
    for lvl in TARGET_LEVELS:
        n = (meta[lvl].notna() & (meta[lvl] != '')).sum()
        print(f"  {lvl}: {n}/{len(meta)}")


if __name__ == '__main__':
    main()
