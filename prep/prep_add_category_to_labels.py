#!/usr/bin/env python3
"""
为 sequence_labels.tsv 追加 Category 列, 并修复 type 列 (细分负样本来源).
逻辑:
  - 正样本: Category = ICTV species (从 selected_viruses.tsv 查找)
  - 负样本: Category = label (即 negative_A/B/C)
  - type 修复: 负样本 type = label (保留 A/B/C 子类型), 正样本 type = 'positive'
"""

import argparse, os, sys, re
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='追加 Category 列 + 修复 type 列')
    parser.add_argument('--labels', required=True, help='sequence_labels.tsv')
    parser.add_argument('--virus-meta', required=True, help='selected_viruses.tsv (含 species 列)')
    parser.add_argument('--out', required=True, help='输出路径')
    args = parser.parse_args()

    labels = pd.read_csv(args.labels, sep='\t')
    meta = pd.read_csv(args.virus_meta, sep='\t')

    # accession -> species 映射
    acc2species = dict(zip(meta['accession'], meta['species']))
    print(f"[meta] {len(acc2species)} accessions mapped")

    # 修复 type 列: 负样本 type = label (保留 negative_A/B/C)
    old_types = labels['type'].value_counts().to_dict()
    labels.loc[labels['label'] != 'positive', 'type'] = labels.loc[labels['label'] != 'positive', 'label']
    new_types = labels['type'].value_counts().to_dict()
    print(f"[type] 修复前: {old_types}")
    print(f"[type] 修复后: {new_types}")

    categories = []
    for _, row in labels.iterrows():
        sid = row['seq_id']
        label = row['label']

        if label == 'positive':
            m = re.search(r'source=([^|]+)', sid)
            acc = m.group(1) if m else None
            if acc and acc in acc2species:
                categories.append(acc2species[acc])
            elif acc:
                categories.append(acc)
            else:
                categories.append(sid[:40])
        else:
            # 负样本: Category = label (negative_A / negative_B / negative_C)
            categories.append(label)

    labels.insert(3, 'Category', categories)

    pos_cats = labels[labels['label'] == 'positive']['Category']
    n_unique = pos_cats.nunique()
    print(f"[category] 正样本 Category 种类: {n_unique}")

    cat_counts = pos_cats.value_counts()
    multi = cat_counts[cat_counts > 1]
    if len(multi) > 0:
        print(f"[category] 多片段实体 (Category 下 >1 条序列): {len(multi)} 个")
        for cat, cnt in multi.sort_values(ascending=False).head(10).items():
            print(f"  {cat}: {cnt} 条序列")

    labels.to_csv(args.out, sep='\t', index=False)
    print(f"[done] {args.out} ({len(labels)} 行, {len(labels.columns)} 列)")

if __name__ == '__main__':
    main()
