#!/usr/bin/env python3
"""
标准化病毒分类预测结果 — 将 MMseqs2/ACVirus/VITAP 统一为:
contig_id | Realm | Kingdom | Phylum | Class | Order | Family | Genus | Species

用法:
  python standardize_classification.py -t mmseqs   -i test_all_lca.tsv                    -o std_mmseqs.tsv
  python standardize_classification.py -t acvirus  -i final_result.tsv                    -o std_acvirus.tsv
  python standardize_classification.py -t vitap    -i best_determined_lineages.tsv        -o std_vitap.tsv
"""

import argparse, csv, os, re, sys

TAX_LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# ── MMseqs2 ──────────────────────────────────────────────────

def parse_lineage_mmseqs(lineage_str):
    """
    解析 MMseqs2 lineage: k_Orthornavirae;p_Duplornaviricota;c_...;s_Fijivirus cuartoense
    映射: -_ → Realm, k_ → Kingdom, p_ → Phylum, c_ → Class, o_ → Order, f_ → Family, g_ → Genus, s_ → Species
    """
    result = {lvl: '' for lvl in TAX_LEVELS}
    if not lineage_str:
        return result
    # 反向匹配: 找 key:value pairs
    # 从 lineage 末尾开始解析，避免 s_ 匹配到 s_ in other names
    rank_map = {'-': 'Realm', 'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class',
                'o': 'Order', 'f': 'Family', 'g': 'Genus', 's': 'Species'}
    # 分割: 以 ; 但注意 s_Fijivirus cuartoense 中间可能有空格
    # 模式: [kpos_-]key_value → 空格是 ;
    parts = []
    for segment in re.split(r';(?=[kpcogfs-]_)', lineage_str):
        segment = segment.strip()
        if '_' not in segment:
            continue
        rank_letter = segment[0]
        if rank_letter in rank_map:
            value = segment[2:] if segment[1] == '_' else segment[2:]  # skip k_ or -_
            rank_name = rank_map[rank_letter]
            if value and value != 'NA':
                result[rank_name] = value
    return result


def standardize_mmseqs(input_path, output_path):
    """MMseqs2 LCA TSV — 无header，col1=seq_id, last_col=lineage"""
    rows = []
    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            contig_id = parts[0].strip()
            lineage_str = parts[-1]  # 最后列
            taxa = parse_lineage_mmseqs(lineage_str)
            rows.append([contig_id] + [taxa.get(l, '') for l in TAX_LEVELS])

    write_output(rows, output_path)
    return len(rows)


# ── ACVirus ──────────────────────────────────────────────────

def standardize_acvirus(input_path, output_path):
    """ACVirus final_result.tsv — header: Nucleotide,Realm,...,Suborder,Family,Subfamily,Genus,Subgenus,Species"""
    rows = []
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig_id = row.get('Nucleotide', '').strip()
            if not contig_id:
                continue
            taxa = []
            for lvl in TAX_LEVELS:
                v = row.get(lvl, '').strip()
                if v in ('', 'NA', 'none', 'Unassigned'):
                    v = ''
                taxa.append(v)
            rows.append([contig_id] + taxa)

    write_output(rows, output_path)
    return len(rows)


# ── VITAP ────────────────────────────────────────────────────

def parse_lineage_vitap(lineage_str):
    """
    解析 VITAP lineage: Riboviria, Orthornavirae, Duplornaviricota, Resentoviricetes, Reovirales, Spinareoviridae, Fijivirus, Fijivirus cuartoense
    逗号分隔，顺序: Realm→Species
    """
    result = {lvl: '' for lvl in TAX_LEVELS}
    if not lineage_str:
        return result
    parts = [p.strip() for p in lineage_str.split(';')]
    for i, lvl in enumerate(TAX_LEVELS):
        if i < len(parts) and parts[i] and parts[i] not in ('NA', '', 'none', 'Unassigned'):
            result[lvl] = parts[i]
    return result


def standardize_vitap(input_path, output_path):
    """VITAP best_determined_lineages.tsv — Genome_ID, lineage, ..."""
    rows = []
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig_id = row.get('Genome_ID', row.get('genome_id', '')).strip()
            if not contig_id:
                continue
            lineage_str = row.get('lineage', row.get('Lineage', '')).strip()
            taxa = parse_lineage_vitap(lineage_str)
            rows.append([contig_id] + [taxa.get(l, '') for l in TAX_LEVELS])

    write_output(rows, output_path)
    return len(rows)


# ── 输出 ─────────────────────────────────────────────────────

def write_output(rows, output_path):
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['contig_id'] + TAX_LEVELS)
        w.writerows(rows)


# ── CLI ──────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='标准化病毒分类结果')
    parser.add_argument('-t', '--tool', required=True, choices=['mmseqs', 'acvirus', 'vitap'])
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[ERROR] File not found: {args.input}")
        sys.exit(1)

    if args.tool == 'mmseqs':
        n = standardize_mmseqs(args.input, args.output)
    elif args.tool == 'acvirus':
        n = standardize_acvirus(args.input, args.output)
    elif args.tool == 'vitap':
        n = standardize_vitap(args.input, args.output)

    print(f"[{args.tool}] Standardized {n} records → {args.output}")


if __name__ == '__main__':
    main()
