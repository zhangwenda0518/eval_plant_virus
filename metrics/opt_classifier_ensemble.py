#!/usr/bin/env python3
"""
分类工具组合最优判断 — 层级共识 + 可靠度门禁 + 文件名/引号容错

用法:
  python opt_classifier_ensemble.py \
      --predictions step9_classification/integrated/ \
      --meta step4_classification_eval/test_metadata_full.tsv \
      --weight-from step9_classification/analysis/ACVirus/stratified_family.tsv \
      --min-tool-accuracy 50.0 \
      --outdir step9_classification/analysis/ensemble_opt/
"""

import argparse, os, sys, itertools, glob
from collections import defaultdict, Counter
import pandas as pd
import numpy as np

TAX_LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
TOOL_LABELS = {
    'acvirus': 'ACVirus', 'mmseqs': 'MMseqs2',
    'vitap': 'VITAP', 'vcontact3': 'vConTACT3',
}
PRIORITY_ORDER = ['vcontact3', 'acvirus', 'vitap', 'mmseqs']
DECOY_TYPES = {"pfam", "eve", "host"}


def _clean(s):
    """脱去外侧双引号和空白"""
    return str(s).strip('"').strip() if s is not None and not pd.isna(s) else ''


def load_tool_predictions(pred_dir, tool_names):
    """文件名大小写自适应 + 双引号净化"""
    file_map = {}
    for f in glob.glob(os.path.join(pred_dir, 'standardized_*.tsv')):
        key = os.path.basename(f).replace('standardized_', '').replace('.tsv', '').lower()
        file_map[key] = f

    tool_data = {}
    for t in tool_names:
        key = t.lower()
        path = file_map.get(key)
        if not path:
            print(f"[WARN] {t}: file not found in {pred_dir}, skipping")
            continue

        df = pd.read_csv(path, sep='\t')
        preds = {}
        for _, row in df.iterrows():
            sid = _clean(row.get('contig_id') or row.get('seq_id') or row.get('id'))
            if not sid:
                continue
            lineage = []
            for lvl in TAX_LEVELS:
                val = _clean(row.get(lvl, row.get(lvl.lower(), row.get(lvl.upper(), ''))))
                if val and val.lower() not in ('na', 'none', 'unassigned', 'nan', ''):
                    lineage.append(val)
                else:
                    lineage.append(None)
            preds[sid] = lineage

        label = TOOL_LABELS.get(key, t)
        tool_data[label] = preds
        print(f"[load] {label}: {len(preds)} contigs ({os.path.basename(path)})")
    return tool_data


def load_ground_truth(meta_path):
    """真值：ID 和分类名均脱引号"""
    meta = pd.read_csv(meta_path, sep='\t')
    truth, types = {}, {}
    for _, row in meta.iterrows():
        sid = _clean(row.get('seq_id') or row.get('source_accession') or row.get('id'))
        if not sid:
            continue
        lineage = []
        for lvl in TAX_LEVELS:
            val = _clean(row.get(lvl, row.get(lvl.lower(), '')))
            if val and val.lower() not in ('na', 'none', 'unassigned', 'nan', ''):
                lineage.append(val)
            else:
                lineage.append(None)
        truth[sid] = lineage
        types[sid] = _clean(row.get('virus_type', 'unknown'))
    print(f"[truth] {len(truth)} sequences")
    return truth, types


def load_stratified_accuracy(path):
    """从目录或单文件读取各工具各等级 OVERALL 准确率"""
    if not path or not os.path.exists(path):
        return None
    try:
        # 支持目录: 扫描所有 stratified_*.tsv
        if os.path.isdir(path):
            files = glob.glob(os.path.join(path, 'stratified_*.tsv'))
            if not files:
                print(f"[WARN] No stratified_*.tsv found in {path}")
                return None
            dfs = []
            for f in files:
                dfs.append(pd.read_csv(f, sep='\t'))
            df = pd.concat(dfs)
        else:
            df = pd.read_csv(path, sep='\t')

        cols = {c.lower(): c for c in df.columns}
        if 'level' not in cols or 'method' not in cols or 'accuracy' not in cols:
            return None
        lvl_col, met_col, acc_col = cols['level'], cols['method'], cols['accuracy']
        stratum_col = cols.get('stratum')
        w = {}
        for _, row in df.iterrows():
            if stratum_col and str(row[stratum_col]).strip().upper() != 'OVERALL':
                continue
            lvl = str(row[lvl_col]).strip().capitalize()
            if lvl not in TAX_LEVELS:
                continue
            raw = str(row[met_col]).strip().lower()
            label = None
            for k, v in TOOL_LABELS.items():
                if k.lower() == raw or v.lower() == raw:
                    label = v; break
            if not label:
                label = row[met_col]
            v = row[acc_col]
            if pd.notna(v):
                w.setdefault(lvl, {})[label] = float(v)
        print(f"[weights] Loaded per-level accuracy for {len(w)} levels "
              f"({', '.join(sorted(w.keys()))}) from {path}")
        return w
    except Exception as e:
        print(f"[WARN] Failed to load weights from {path}: {e}")
        return None


def get_hierarchical_consensus(tool_lineages, strategy, subset_tools,
                                level_weights=None, min_tool_acc=0.0):
    """自上而下层级共识 + 不可靠工具单兵深入门禁"""
    path = []
    n_total = len(subset_tools)

    for i, lvl_name in enumerate(TAX_LEVELS):
        compat = [t for t in subset_tools
                  if t in tool_lineages and tool_lineages[t][:i] == path[:i]]
        if not compat:
            break

        votes = {t: tool_lineages[t][i] for t in compat
                 if tool_lineages[t][i] is not None}
        if not votes:
            break

        vals = list(votes.values())
        counter = Counter(vals)
        top = counter.most_common(2)
        top_val, top_cnt = top[0]
        has_conflict = len(counter) > 1
        accepted = None

        if strategy == 'any_priority':
            for t in subset_tools:
                if t in votes:
                    accepted = votes[t]; break
        elif strategy == 'unanimous':
            if len(votes) == n_total and not has_conflict:
                accepted = top_val
        elif strategy == 'majority':
            if has_conflict:
                if top_cnt >= n_total / 2.0 and (len(top) == 1 or top_cnt > top[1][1]):
                    accepted = top_val
            else:
                accepted = top_val
        elif strategy == 'majority2':
            if has_conflict:
                if top_cnt >= 2 and (len(top) == 1 or top_cnt > top[1][1]):
                    accepted = top_val
            else:
                accepted = top_val
        elif strategy == 'weighted' and level_weights:
            w = level_weights.get(lvl_name, {})
            vote_w = defaultdict(float)
            for t, v in votes.items():
                vote_w[v] += w.get(t, 1.0)
            top_w = sorted(vote_w.items(), key=lambda x: x[1], reverse=True)
            if len(top_w) > 1 and abs(top_w[0][1] - top_w[1][1]) < 0.001:
                accepted = None
            else:
                accepted = top_w[0][0]

        if not accepted:
            break

        if not has_conflict and len(votes) == 1:
            proposer = list(votes.keys())[0]
            if level_weights and lvl_name in level_weights:
                hist_acc = level_weights[lvl_name].get(proposer, 100)
                if hist_acc < min_tool_acc:
                    break

        path.append(accepted)

    while len(path) < len(TAX_LEVELS):
        path.append(None)
    return path


def evaluate_ensemble(truth, types, tool_data, subset_tools, strategy,
                       level_weights=None, min_tool_acc=0.0):
    preds = {}
    for sid in truth:
        tl = {t: tool_data[t][sid] for t in subset_tools
              if t in tool_data and sid in tool_data[t]}
        preds[sid] = get_hierarchical_consensus(tl, strategy, subset_tools,
                                                 level_weights, min_tool_acc)

    metrics = {}
    for lvl_idx, lvl in enumerate(TAX_LEVELS):
        known = dict(total=0, assigned=0, correct=0)
        novel = dict(total=0, assigned=0, correct=0)
        decoy = dict(total=0, assigned=0)
        for sid, true_lin in truth.items():
            vt = types.get(sid, '')
            pv = preds[sid][lvl_idx]
            if vt in DECOY_TYPES:
                decoy['total'] += 1
                if pv is not None:
                    decoy['assigned'] += 1
            else:
                tv = true_lin[lvl_idx]
                if tv is None:
                    continue
                d = known if vt == 'known' else novel
                d['total'] += 1
                if pv is not None:
                    d['assigned'] += 1
                    if pv == tv:
                        d['correct'] += 1

        out = {}
        for label, d in [('known', known), ('novel', novel)]:
            t, a, c = d['total'], d['assigned'], d['correct']
            out[label] = (t, a, c,
                          (c / a * 100) if a > 0 else None,
                          (a / t * 100) if t > 0 else 0.0)
        out['decoy'] = (decoy['total'], decoy['assigned'],
                        decoy['assigned'] / decoy['total'] * 100 if decoy['total'] > 0 else 0.0)
        metrics[lvl] = out
    return metrics


def sort_by_priority(labels):
    ll = {lbl.lower(): lbl for lbl in labels}
    order = [ll[p] for p in PRIORITY_ORDER if p in ll]
    order += [lbl for lbl in labels if lbl not in order]
    return order


def main():
    parser = argparse.ArgumentParser(description='最优组合（层级共识+门禁+净化版）')
    parser.add_argument('--predictions', default='.')
    parser.add_argument('--meta', required=True)
    parser.add_argument('--tools', default=None)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--min-assign', type=float, default=5.0)
    parser.add_argument('--weight-from', default=None,
                        help='stratified_*.tsv，提取各工具各等级历史准确率')
    parser.add_argument('--min-tool-accuracy', type=float, default=50.0,
                        help='无冲突单工具深入时的历史准确率最低阈值(%)')
    parser.add_argument('--sort-by', default='combined_score',
                        choices=['combined_score', 'f1_score', 'accuracy'],
                        help='排序指标: combined_score (默认,F1×抗FP), f1_score, accuracy')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.tools:
        tool_names = [t.strip() for t in args.tools.split(',')]
    else:
        tool_names = sorted(set(
            os.path.basename(f).replace('standardized_', '').replace('.tsv', '')
            for f in glob.glob(os.path.join(args.predictions, 'standardized_*.tsv'))))
    print(f"Tools ({len(tool_names)}): {', '.join(tool_names)}")
    tool_data = load_tool_predictions(args.predictions, tool_names)
    truth, types = load_ground_truth(args.meta)

    level_weights = load_stratified_accuracy(args.weight_from)
    if level_weights and args.min_tool_accuracy > 0:
        print(f"[gate] Reject uncontradicted deepening if tool accuracy < {args.min_tool_accuracy}%")
    elif args.min_tool_accuracy > 0:
        print("[WARN] --min-tool-accuracy is set but no valid weights file loaded. Reliability gate is inactive.")

    all_labels = list(tool_data.keys())
    strategies = ['unanimous', 'majority', 'majority2', 'any_priority']
    if level_weights:
        strategies.append('weighted')

    def make_row(subset, strategy, lvl, out, vt, n):
        t, a, c, acc, ar = out[vt]
        fpr = out['decoy'][2]
        p = (c / a * 100) if a > 0 else 0
        r = (c / t * 100) if t > 0 else 0
        f1 = round(2 * p * r / (p + r), 1) if (p + r) > 0 else None
        cs = round(f1 * max(0, 1 - fpr / 100), 1) if f1 else None
        return {
            'subset': subset, 'strategy': strategy, 'level': lvl,
            'virus_type': vt, 'total': t, 'assigned': a, 'correct': c,
            'accuracy': round(acc, 1) if acc else None,
            'assignment_rate': round(ar, 1),
            'f1_score': f1, 'fp_rate': round(fpr, 1),
            'fp_assigned': out['decoy'][1], 'fp_total': out['decoy'][0],
            'combined_score': cs, 'n_tools': n,
        }

    all_rows = []
    for lbl in all_labels:
        res = evaluate_ensemble(truth, types, tool_data, [lbl], 'any_priority',
                                 level_weights, args.min_tool_accuracy)
        for lvl in TAX_LEVELS:
            for vt in ['known', 'novel']:
                all_rows.append(make_row(lbl, 'single', lvl, res[lvl], vt, 1))

    total = sum(1 for n in range(2, len(all_labels) + 1)
                for _ in itertools.combinations(all_labels, n))
    done = 0
    for n in range(2, len(all_labels) + 1):
        for subset in itertools.combinations(all_labels, n):
            ordered = sort_by_priority(list(subset))
            name = '+'.join(ordered)
            done += 1
            if done % 20 == 0 or done == total:
                print(f"  [{done}/{total}] {name} ...")
            valid = [s for s in strategies
                     if not (s == 'majority' and n == 2)
                     and not (s == 'majority2' and n < 3)]
            for s in valid:
                res = evaluate_ensemble(truth, types, tool_data, ordered, s,
                                         level_weights, args.min_tool_accuracy)
                for lvl in TAX_LEVELS:
                    for vt in ['known', 'novel']:
                        all_rows.append(make_row(name, s, lvl, res[lvl], vt, n))

    df = pd.DataFrame(all_rows)
    df = df[df['assignment_rate'] >= args.min_assign]
    df.to_csv(os.path.join(args.outdir, 'ensemble_optimization.tsv'), sep='\t', index=False, na_rep='NA')

    sc = args.sort_by
    print(f"\n{'='*90}")
    print(f"  最优组合（层级共识 + 门禁={args.min_tool_accuracy}% + 排序={sc}）")

    summary_rows = []
    for vt_label, vt_title in [('known', '已知病毒'), ('novel', '新病毒')]:
        vdf = df[df['virus_type'] == vt_label]
        print(f"\n  ── {vt_title} ──")
        print(f"  {'Level':10s} | {'Best subset':28s} | {'Strategy':12s} | {'Acc':>5s} | {'Assign':>6s} | {'F1':>4s} | {'FP%':>5s} | {'CS':>4s} | N")
        print(f"  {'-'*10}-+-{'-'*28}-+-{'-'*12}-+-{'-'*5}-+-{'-'*6}-+-{'-'*4}-+-{'-'*5}-+-{'-'*4}-+--")
        for lvl in TAX_LEVELS:
            ldf = vdf[vdf['level'] == lvl].sort_values(sc, ascending=False, na_position='last')
            if ldf.empty: continue
            b = ldf.iloc[0]
            summary_rows.append(b.to_dict())

            acc_str = f"{b['accuracy']:.0f}%" if pd.notna(b['accuracy']) else "NA"
            f1_str = f"{b['f1_score']:.0f}" if pd.notna(b['f1_score']) else "NA"
            fp_str = f"{b['fp_rate']:.0f}%" if pd.notna(b['fp_rate']) else "NA"
            cs_str = f"{b['combined_score']:.0f}" if pd.notna(b['combined_score']) else "NA"
            ar_str = f"{b['assignment_rate']:.1f}%" if pd.notna(b['assignment_rate']) else "NA"

            print(f"  {lvl:10s} | {b['subset']:28s} | {b['strategy']:12s} | {acc_str:>5} | {ar_str:>6} | {f1_str:>4} | {fp_str:>5} | {cs_str:>4} | {b['n_tools']:>2d}")

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(args.outdir, 'ensemble_best_per_level.tsv')
        summary_df.to_csv(summary_path, sep='\t', index=False, na_rep='NA')
        print(f"\n[save] Best configurations per level saved to: {summary_path}")

    print(f"\nDone. {args.outdir}")


if __name__ == '__main__':
    main()
