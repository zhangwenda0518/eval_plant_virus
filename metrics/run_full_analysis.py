#!/usr/bin/env python3
"""
分类评估全面分析 — 一键运行所有维度的分层评估并生成汇总报告

维度:
  1. 突变率影响 (stratum=mut_type)                  → 各工具的突变耐受度
  2. 覆盖度影响 (stratum=cov_type)                  → 各工具的覆盖度鲁棒性
  3. 病毒类型对比 (stratum=virus_type)              → known vs novel vs decoy
  4. 突变×覆盖交互 (filter mut_type=stratum cov_type) → 最差场景叠加
  5. decoy 假阳性 (filter virus_type=pfam,eve,host) → 各工具的假阳性率
  6. 突变×病毒类型 (filter virus_type=stratum mut_type) → 新病毒突变脆弱性

用法:
  python run_full_analysis.py \
      --predictions-dir step9_classification/integrated/ \
      --meta step4_classification_eval/test_metadata_full.tsv \
      --tools acvirus,mmseqs,vitap,vcontact3 \
      --outdir step9_classification/analysis/
"""

import argparse, os, sys, subprocess


def run_cmd(cmd):
    print(f"\n{'='*60}")
    print(f"  RUN: {cmd}")
    print(f"{'='*60}")
    ret = os.system(cmd)
    if ret != 0:
        print(f"[WARN] Exit code {ret}: {cmd}")
    return ret


def main():
    parser = argparse.ArgumentParser(description='分类评估全面分析')
    parser.add_argument('--predictions-dir', required=True,
                        help='标准化预测文件目录')
    parser.add_argument('--meta', required=True)
    parser.add_argument('--tools', default=None,
                        help='指定工具（逗号分隔），不指定则自动扫描 standardized_*.tsv')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    script = 'eval_plant_virus/metrics/calc_classification_stratified.py'

    if args.tools:
        all_tools = [t.strip() for t in args.tools.split(',')]
    else:
        import glob
        pred_dir = args.predictions_dir
        all_tools = sorted(set(
            os.path.basename(f).replace('standardized_', '').replace('.tsv', '')
            for f in glob.glob(os.path.join(pred_dir, 'standardized_*.tsv'))
        ))
    if not all_tools:
        print(f"[ERROR] No tools found in {args.predictions_dir}")
        sys.exit(1)
    print(f"Tools ({len(all_tools)}): {', '.join(all_tools)}")

    # ============================================================
    # Phase 1: 各工具全部维度（不筛选，auto-strata）
    # ============================================================
    print("\n" + "=" * 70)
    print("  Phase 1: 各工具全部维度 (mut_type + cov_type + virus_type)")
    print("=" * 70)

    for tool in all_tools:
        pred_file = os.path.join(args.predictions_dir, f'standardized_{tool}.tsv')
        if not os.path.exists(pred_file):
            print(f"  [SKIP] {tool}: prediction file not found ({pred_file})")
            continue
        out_dir = os.path.join(args.outdir, tool)
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--outdir {out_dir}')

    # ============================================================
    # Phase 2: 维度专项分析
    # ============================================================
    print("\n" + "=" * 70)
    print("  Phase 2: 维度专项分析")
    print("=" * 70)

    for tool in all_tools:
        pred_file = os.path.join(args.predictions_dir, f'standardized_{tool}.tsv')
        if not os.path.exists(pred_file):
            continue

        # 2a. 已知病毒的突变率影响
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=known --stratum mut_type '
                f'--outdir {args.outdir}/{tool}/known_by_mut/')

        # 2b. 已知病毒的覆盖度影响
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=known --stratum cov_type '
                f'--outdir {args.outdir}/{tool}/known_by_cov/')

        # 2c. 新病毒的突变率影响
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=novel --stratum mut_type '
                f'--outdir {args.outdir}/{tool}/novel_by_mut/')

        # 2d. 新病毒的覆盖度影响
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=novel --stratum cov_type '
                f'--outdir {args.outdir}/{tool}/novel_by_cov/')

        # 2e. decoy 假阳性分析
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=pfam,eve,host '
                f'--outdir {args.outdir}/{tool}/decoys/')

    # ============================================================
    # Phase 3: 最差场景 — 高突变 × 低覆盖叠加
    # ============================================================
    print("\n" + "=" * 70)
    print("  Phase 3: 最差场景 (高突变 × 低覆盖)")
    print("=" * 70)

    for tool in all_tools:
        pred_file = os.path.join(args.predictions_dir, f'standardized_{tool}.tsv')
        if not os.path.exists(pred_file):
            continue

        # 3a. mut15 时覆盖度影响
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=known --filter mut_type=mut15 '
                f'--stratum cov_type '
                f'--outdir {args.outdir}/{tool}/worst_mut15_by_cov/')

        # 3b. cov40 时突变率影响
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=known --filter cov_type=cov40 '
                f'--stratum mut_type '
                f'--outdir {args.outdir}/{tool}/worst_cov40_by_mut/')

        # 3c. novel + mut15 + 覆盖度
        run_cmd(f'python {script} --predictions {pred_file} --meta {args.meta} '
                f'--filter virus_type=novel --filter mut_type=mut15 '
                f'--stratum cov_type '
                f'--outdir {args.outdir}/{tool}/worst_novel_mut15/')

    # ============================================================
    # Phase 4: 汇总可视化
    # ============================================================
    print("\n" + "=" * 70)
    print("  Phase 4: 汇总可视化")
    print("=" * 70)

    plot_script = 'eval_plant_virus/metrics/plot_classification_comparison.py'
    run_cmd(f'python {plot_script} --input {args.outdir} '
            f'--outdir {args.outdir}/plots/')

    # 子维度可视化
    sub_dims = ['known_by_mut', 'known_by_cov', 'novel_by_mut', 'novel_by_cov', 'decoys']
    for dim in sub_dims:
        dim_dir = os.path.join(args.outdir, dim)
        if os.path.isdir(dim_dir):
            run_cmd(f'python {plot_script} --input {dim_dir} '
                    f'--outdir {args.outdir}/plots/{dim}/')

    # ============================================================
    # Phase 5: 最优工具组合
    # ============================================================
    print("\n" + "=" * 70)
    print("  Phase 5: 最优工具组合")
    print("=" * 70)

    opt_script = 'eval_plant_virus/metrics/opt_classifier_ensemble.py'
    run_cmd(f'python {opt_script} --predictions {args.predictions_dir} '
            f'--meta {args.meta} '
            f'--outdir {args.outdir}/ensemble_opt/ '
            f'--filter-virus-type known,novel')

    print(f"\n{'='*70}")
    print(f"  ALL DONE")
    print(f"{'='*70}")
    print(f"  Output: {args.outdir}")
    print(f"  Plots:  {args.outdir}/plots/")
    print(f"  Summary text: {args.outdir}/plots/summary_all_levels.txt")


if __name__ == '__main__':
    main()
