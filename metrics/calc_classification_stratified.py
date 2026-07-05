#!/usr/bin/env python3
"""
病毒分类评估分层报告 v3 — 支持多维筛选与分层

分层维度: virus_type | mut_type | cov_type | coverage_pct | mutation_rate_pct
筛选: --filter virus_type=novel --filter mut_type=mut5

用法:
  # 按 virus_type 分层
  python calc_classification_stratified.py --stratum virus_type ...

  # 只看 mut5 的 known 病毒
  python calc_classification_stratified.py \
      --filter virus_type=known --filter mut_type=mut5 --stratum cov_type ...

  # 只看 decoy (pfam/eve/host)，计算假阳性率
  python calc_classification_stratified.py \
      --filter virus_type=pfam --filter virus_type=eve --filter virus_type=host ...
"""

import argparse, os
import pandas as pd
import numpy as np
from collections import defaultdict

DECOY_TYPES = {"pfam", "eve", "host"}


def load_predictions(pred_tsv):
    return pd.read_csv(pred_tsv, sep='\t')


def load_metadata(meta_tsv):
    return pd.read_csv(meta_tsv, sep='\t')


def detect_stratum(meta_df, specified=None):
    """解析分层列（支持逗号分隔的多列）"""
    if specified and specified in meta_df.columns:
        return [specified]
    if specified and ',' in specified:
        return [c.strip() for c in specified.split(',') if c.strip() in meta_df.columns]
    # 不指定则默认三维交叉
    return [c for c in ['virus_type', 'mut_type', 'cov_type'] if c in meta_df.columns]


def stratum_label(meta_row, stratum_cols):
    """返回分层标签（多列用 | 连接）"""
    if isinstance(stratum_cols, str):
        stratum_cols = [stratum_cols]
    parts = []
    for col in stratum_cols:
        raw = meta_row.get(col, '')
        if raw is None or str(raw).strip() == '' or pd.isna(raw):
            parts.append('decoy')
        elif col in ('mutation_rate_pct', 'coverage_pct'):
            parts.append(f"{int(float(raw))}%")
        else:
            parts.append(str(raw).strip())
    return '|'.join(parts)


def stratum_sort_key(s):
    try:
        return -float(s.replace('%', '').replace('|', ''))
    except (ValueError, AttributeError):
        return 0


def compute_stratified_metrics_fast(pred_df, meta_df, level, stratum_cols, pred_map):
    """快速版：pred_map 和 stratum_labels 已预构建"""
    stratum_col_str = ','.join(stratum_cols)
    # 预取 strata 和 true_vals
    meta_sids = meta_df['_sid'].values
    meta_strata = meta_df['_stratum_' + stratum_col_str].values
    meta_true = meta_df['_true_' + level].values
    meta_is_decoy = meta_df['_is_decoy'].values

    results = defaultdict(lambda: {"correct": 0, "total": 0, "unassigned": 0, "decoy_assigned": 0})

    for i in range(len(meta_sids)):
        sid = meta_sids[i]
        stratum = meta_strata[i]
        true_val = meta_true[i]
        is_decoy = meta_is_decoy[i]

        pred_val = pred_map.get(sid, None)

        results[stratum]["total"] += 1
        if pred_val is None:
            results[stratum]["unassigned"] += 1
        elif is_decoy:
            results[stratum]["decoy_assigned"] += 1
        elif pred_val == true_val:
            results[stratum]["correct"] += 1

    print(f"\n{'='*60}")
    print(f"  Level: {level.upper()}  |  Stratum: {stratum_col_str}")
    print(f"{'='*60}")
    has_decoy = any(r.get("decoy_assigned", 0) > 0 for r in results.values())
    if has_decoy:
        print(f"  {'Stratum':12s} | {'Assign':>9s} | {'Acc':>6s} | {'Error':>6s} | {'Miss':>6s} | {'FP':>6s}")
        print(f"  {'-'*12}-+-{'-'*9}-+-{'-'*6}-+-{'-'*6}-+-{'-'*6}-+-{'-'*6}")
    else:
        print(f"  {'Stratum':12s} | {'Assign':>9s} | {'Acc':>6s} | {'Error':>6s} | {'Miss':>6s}")
        print(f"  {'-'*12}-+-{'-'*9}-+-{'-'*6}-+-{'-'*6}-+-{'-'*6}")

    rows = []
    strata_order = sorted(results.keys(), key=stratum_sort_key)
    for stratum in strata_order:
        r = results[stratum]
        if r['total'] == 0:
            continue
        assigned = r['total'] - r['unassigned']
        effective_correct = r['correct']
        acc = effective_correct / assigned * 100 if assigned > 0 else None
        err = (assigned - effective_correct) / assigned * 100 if assigned > 0 else None
        assign_rate = assigned / r['total'] * 100 if r['total'] > 0 else 0
        miss_rate = r['unassigned'] / r['total'] * 100 if r['total'] > 0 else 0
        fp_rate = (r.get('decoy_assigned', 0) / assigned * 100 if assigned > 0 else 0) if has_decoy else None

        acc_str = f"{acc:4.1f}%" if acc is not None else "  N/A"
        err_str = f"{err:4.1f}%" if err is not None else "  N/A"
        if has_decoy:
            fp_str = f"{fp_rate:4.1f}%" if fp_rate is not None else "  N/A"
            print(f"  {stratum:12s} | {assign_rate:4.1f}% ({assigned:3d}/{r['total']:3d}) "
                  f"| {acc_str} | {err_str} "
                  f"| {miss_rate:4.1f}% ({r['unassigned']:3d})"
                  f" | {fp_str}")
        else:
            print(f"  {stratum:12s} | {assign_rate:4.1f}% ({assigned:3d}/{r['total']:3d}) "
                  f"| {acc_str} | {err_str} "
                  f"| {miss_rate:4.1f}% ({r['unassigned']:3d})")

        row = {
            'level': level, 'stratum': stratum, 'stratum_col_str': stratum_col_str,
            'total': r['total'], 'assigned': assigned,
            'unassigned': r['unassigned'], 'correct': effective_correct,
            'accuracy': round(acc, 2) if acc is not None else None,
            'error_rate': round(err, 2) if err is not None else None,
            'assignment_rate': round(assign_rate, 2),
            'missing_rate': round(miss_rate, 2),
        }
        if has_decoy:
            row['decoy_assigned'] = r.get('decoy_assigned', 0)
            row['fp_rate'] = round(fp_rate, 2) if fp_rate is not None else None
        rows.append(row)

    # Overall
    total_c = sum(r['correct'] for r in results.values())
    total_a = sum(r['total'] for r in results.values())
    total_u = sum(r['unassigned'] for r in results.values())
    total_decoy_fp = sum(r.get('decoy_assigned', 0) for r in results.values())
    assigned = total_a - total_u
    overall_acc = total_c / assigned * 100 if assigned > 0 else None
    overall_err = (assigned - total_c) / assigned * 100 if assigned > 0 else None
    overall_assign = assigned / total_a * 100 if total_a > 0 else 0
    overall_miss = total_u / total_a * 100 if total_a > 0 else 0
    overall_fp = total_decoy_fp / assigned * 100 if assigned > 0 else 0

    acc_str = f"{overall_acc:4.1f}%" if overall_acc is not None else "  N/A"
    err_str = f"{overall_err:4.1f}%" if overall_err is not None else "  N/A"
    if has_decoy:
        print(f"  {'OVERALL':12s} | {overall_assign:4.1f}% ({assigned:3d}/{total_a:3d}) "
              f"| {acc_str} | {err_str} "
              f"| {overall_miss:4.1f}% ({total_u:3d})"
              f" | {overall_fp:4.1f}%")
    else:
        print(f"  {'OVERALL':12s} | {overall_assign:4.1f}% ({assigned:3d}/{total_a:3d}) "
              f"| {acc_str} | {err_str} "
              f"| {overall_miss:4.1f}% ({total_u:3d})")

    row = {
        'level': level, 'stratum': 'OVERALL', 'stratum_col_str': stratum_col_str,
        'total': total_a, 'assigned': assigned,
        'unassigned': total_u, 'correct': total_c,
        'accuracy': round(overall_acc, 2) if overall_acc is not None else None,
        'error_rate': round(overall_err, 2) if overall_err is not None else None,
        'assignment_rate': round(overall_assign, 2),
        'missing_rate': round(overall_miss, 2),
    }
    if has_decoy:
        row['decoy_assigned'] = total_decoy_fp
        row['fp_rate'] = round(overall_fp, 2)
    rows.append(row)
    return rows


def apply_filters(meta_df, filters):
    """应用多维筛选: --filter virus_type=novel --filter mut_type=mut5"""
    if not filters:
        return meta_df
    for f in filters:
        if '=' in f:
            col, val = f.split('=', 1)
        else:
            # 兼容旧格式: --filter-type X → 等同于 --filter virus_type=X
            col, val = 'virus_type', f
        if col in meta_df.columns:
            if ',' in val:
                meta_df = meta_df[meta_df[col].isin(val.split(','))]
            else:
                meta_df = meta_df[meta_df[col] == val]
            print(f"[filter] {col}={val} → {len(meta_df)} rows")
        else:
            print(f"[filter] WARNING: column '{col}' not found in metadata")
    return meta_df


def main():
    parser = argparse.ArgumentParser(description="病毒分类评估分层报告 v3")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--meta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--stratum", default=None,
                        help="分层列: virus_type / mut_type / cov_type / coverage_pct / mutation_rate_pct")
    parser.add_argument("--filter", action="append", default=[],
                        help="筛选条件 col=val (可重复多次)")
    parser.add_argument("--filter-type", default=None,
                        help="[兼容旧版] 等同 --filter virus_type=VALUE")
    parser.add_argument("--levels", default=None,
                        help="分类等级 (逗号分隔)，默认自动检测")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    pred_df = load_predictions(args.predictions)
    meta_df = load_metadata(args.meta)

    # 筛选
    filters = args.filter.copy()
    if args.filter_type:
        filters.append(args.filter_type)  # 旧兼容
    meta_df = apply_filters(meta_df, filters)
    if len(meta_df) == 0:
        print("[ERROR] No rows after filtering!")
        return

    # 确定分层维度：指定则用指定，否则三维交叉 + 各单维度
    ALL_STRATA = ['mut_type', 'cov_type', 'virus_type',
                  'mutation_rate_pct', 'coverage_pct']
    if args.stratum:
        # 手工指定：逗号分隔=多列交叉，单列=单维度
        cols = detect_stratum(meta_df, args.stratum)
        strata = [cols]  # 一个分组（可能是多列）
    else:
        # 默认：三维交叉 + 各单维度
        cross_cols = [c for c in ['virus_type', 'mut_type', 'cov_type'] if c in meta_df.columns]
        strata = [cross_cols]  # 三维交叉
        # 再加各个单维度
        for c in ALL_STRATA:
            if c in meta_df.columns and c not in cross_cols:
                strata.append([c])
    print(f"Strata: {len(strata)} dimensions  |  Rows: {len(meta_df)}")

    # 检测可用的分类等级
    pred_cols_lower = [c.lower() for c in pred_df.columns]
    meta_cols_lower = [c.lower() for c in meta_df.columns]
    if args.levels:
        level_order = [l.strip().lower() for l in args.levels.split(',')]
    else:
        level_order = ['realm', 'kingdom', 'phylum', 'class', 'order',
                       'family', 'genus', 'species']
    levels = [l for l in level_order if l in pred_cols_lower and l in meta_cols_lower]
    if not levels:
        levels = ['family', 'genus', 'species']
    skipped = [l for l in level_order if l not in levels]
    print(f"Levels: {levels}")
    if skipped:
        print(f"Skipped: {skipped}")

    # ── 预构建加速结构 ──
    print("Pre-building lookup structures...")
    # pred_maps: {level: {sid: pred_val}}
    pred_maps = {lvl: {} for lvl in levels}
    for _, row in pred_df.iterrows():
        sid = row.get('contig_id') or row.get('seq_id') or row.get('id')
        if not sid:
            continue
        for lvl in levels:
            v = row.get(lvl, row.get(lvl.capitalize(), ''))
            if v and str(v).strip() not in ('', 'NA', 'none', 'Unassigned'):
                pred_maps[lvl][sid] = str(v).strip()

    # 预取 meta 列
    meta_df['_sid'] = [row.get('seq_id') or row.get('source_accession') or row.get('id')
                       for _, row in meta_df.iterrows()]
    meta_df['_is_decoy'] = meta_df['virus_type'].isin(DECOY_TYPES).values
    for lvl in levels:
        meta_df['_true_' + lvl] = [str(row.get(lvl, row.get(lvl.capitalize(), ''))).strip()
                                    for _, row in meta_df.iterrows()]
    for stratum_cols in strata:
        key = ','.join(stratum_cols)
        meta_df['_stratum_' + key] = [stratum_label(row, stratum_cols)
                                       for _, row in meta_df.iterrows()]
    print("  Done.")

    # 遍历分层 × 等级
    total_rows = 0
    for level in levels:
        level_rows = []
        for stratum_cols in strata:
            rows = compute_stratified_metrics_fast(
                pred_df, meta_df, level, stratum_cols, pred_maps[level])
            level_rows.extend(rows)
        out_df = pd.DataFrame(level_rows)
        # 拆分为 stratified（交叉表）+ summary（数值汇总）
        cross_col = ','.join(['virus_type', 'mut_type', 'cov_type'])
        cross_df = out_df[out_df['stratum_col_str'] == cross_col]
        p1 = os.path.join(args.outdir, f"stratified_{level}.tsv")
        cross_df.to_csv(p1, sep='\t', index=False, na_rep='NA')
        print(f"  → {p1} ({len(cross_df)} rows)")

        summ_df = out_df[out_df['stratum_col_str'].isin(['mutation_rate_pct', 'coverage_pct'])]
        p2 = os.path.join(args.outdir, f"summary_{level}.tsv")
        summ_df.to_csv(p2, sep='\t', index=False, na_rep='NA')
        print(f"  → {p2} ({len(summ_df)} rows)")

        # 3) 文本格式总结表（screen 输出格式保存为 TSV）
        for label, col_filter in [("all", cross_col),
                                   ("mutation", "mutation_rate_pct"),
                                   ("coverage", "coverage_pct")]:
            sub = out_df[out_df['stratum_col_str'] == col_filter]
            if sub.empty:
                continue
            txt_rows = []
            for _, r in sub.iterrows():
                acc_s = f"{r['accuracy']:.1f}" if not pd.isna(r['accuracy']) else "N/A"
                txt_rows.append({
                    'Stratum': r['stratum'],
                    'Assign_pct': round(r['assignment_rate'], 1),
                    'Assign_n': f"{int(r['assigned'])}/{int(r['total'])}",
                    'Acc_pct': r['accuracy'],
                    'Error_pct': r.get('error_rate', ''),
                    'Miss_pct': round(r['missing_rate'], 1),
                    'Miss_n': int(r['unassigned']),
                    'FP_pct': r.get('fp_rate', ''),
                })
            p3 = os.path.join(args.outdir, f"summary_{level}.Stratum.{label}.tsv")
            pd.DataFrame(txt_rows).to_csv(p3, sep='\t', index=False, na_rep='')
            print(f"  → {p3} ({len(txt_rows)} rows)")
        total_rows += len(cross_df) + len(summ_df)
    print(f"\nSaved {len(levels)} files, {total_rows} total rows → {args.outdir}/")


if __name__ == "__main__":
    main()
