#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Virus Assembly Benchmark — Phase 3: 资源消耗解析
================================================================================

功能: 纯解析，不运行新任务
  自动扫描所有 .time.mem.log 文件，双格式兼容解析，
  输出 resource_summary.tsv

来源:
  1. 基础工具: {sn}.megahit.time.mem.log, {sn}.rnaviralspades.time.mem.log, {sn}.penguin.time.mem.log
  2. 顶层汇总: {sn}.time.mem.log (如果存在)
  3. Merge 工具: {sn}_{subdir}/{sn}_{subdir}.time.mem.log
  4. RefineC_Merge: {sn}_all_tools_refineC_merge/*.time.mem.log

用法:
  python scripts/benchmark_resource.py -i step6_assemblies -o benchmark_results --mode 7

输出:
  - benchmark_results/resource_summary.tsv
"""

import os, sys, argparse, glob

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_utils import (
    discover_samples, parse_resource_log, parse_dir_name,
    BASE_TOOLS, MERGE_TOOL_SUBDIR, REFINEC_SUBDIR_SUFFIX
)


def main():
    parser = argparse.ArgumentParser(description="🧬 Phase 3: 资源消耗解析")
    parser.add_argument("-i", "--input-dir", default="step6_assemblies")
    parser.add_argument("-o", "--out-dir", default="benchmark_results")
    parser.add_argument("--mode", choices=['4', '6', '7', '10'], default='6')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    samples = discover_samples(args.input_dir)
    if not samples:
        print("❌ 未发现任何样本目录！")
        return 1
    print(f"🔍 扫描 {len(samples)} 个样本的资源日志...")

    rows = []
    stats_fmt = {'time-f': 0, 'time-v': 0, 'missing': 0}

    for s in samples:
        sdir = s['sample_dir']
        sn = s['base_name']
        info = {'Mutation_Rate': s.get('Mutation_Rate', '0.0'),
                'Depth': s.get('Depth', 0), 'Rep': s.get('Rep', 'rep1'),
                'Group': s.get('Group', 0),
                'Background_Ratio': s.get('background_ratio', None)}

        # --- 1. 基础组装工具: {sn}.{keyword}.time.mem.{log,csv} ---
        base_keywords = {
            "Megahit":        "megahit",
            "RNAViralSPAdes": "rnaviralspades",
            "Penguin":        "penguin",
            "MetaViralSPAdes":"metaviralspades",
            "Trinity":        "trinity",
            "RNABloom":       "rnabloom",
        }
        for tool, keyword in base_keywords.items():
            found = False
            patterns = [
                os.path.join(sdir, f"{sn}.{keyword}.time.mem.log"),
                os.path.join(sdir, f"*.{keyword}.time.mem.log"),
                os.path.join(sdir, f"{sn}.{keyword}.time.mem.csv"),
                os.path.join(sdir, f"*.{keyword}.time.mem.csv"),
            ]
            for pat in patterns:
                for logfile in glob.glob(pat):
                    res = parse_resource_log(logfile)
                    if res['wall_time_s'] is not None:
                        row = info.copy()
                        row.update({
                            'Tool': tool, 'Wall_Time_s': res['wall_time_s'],
                            'Max_RSS_MB': res['max_rss_mb'], 'CPU_pct': res['cpu_pct'],
                            'IO_In_MB': res['io_in_mb'], 'IO_Out_MB': res['io_out_mb'],
                            'Source_File': os.path.relpath(logfile, sdir),
                            'Source_Format': res['source_format'],
                        })
                        rows.append(row)
                        stats_fmt[res['source_format'] or 'missing'] += 1
                        found = True
                        break
                if found:
                    break

        # --- 2. Merge 工具: {sn}_{subdir}/{sn}_{subdir}.time.mem.log ---
        for tool, subdir in MERGE_TOOL_SUBDIR.items():
            found = False
            patterns = [
                os.path.join(sdir, f"{sn}_{subdir}", f"{sn}_{subdir}.time.mem.log"),
                os.path.join(sdir, f"{sn}_{subdir}", "*.time.mem.log"),
                os.path.join(sdir, f"{sn}_{subdir}", "time.log"),
            ]
            for pat in patterns:
                for logfile in glob.glob(pat):
                    res = parse_resource_log(logfile)
                    if res['wall_time_s'] is not None:
                        row = info.copy()
                        row.update({
                            'Tool': tool, 'Wall_Time_s': res['wall_time_s'],
                            'Max_RSS_MB': res['max_rss_mb'], 'CPU_pct': res['cpu_pct'],
                            'IO_In_MB': res['io_in_mb'], 'IO_Out_MB': res['io_out_mb'],
                            'Source_File': os.path.relpath(logfile, sdir),
                            'Source_Format': res['source_format'],
                        })
                        rows.append(row)
                        stats_fmt[res['source_format'] or 'missing'] += 1
                        found = True
                        break
                if found:
                    break

        # --- 3. RefineC_Merge: {sn}_all_tools_refineC_merge/*.time.mem.log ---
        refinec_dir = os.path.join(sdir, f"{sn}_all_tools_refineC_merge")
        for logfile in glob.glob(os.path.join(refinec_dir, "*.time.mem.log")):
            res = parse_resource_log(logfile)
            if res['wall_time_s'] is not None:
                row = info.copy()
                row.update({
                    'Tool': 'RefineC_Merge', 'Wall_Time_s': res['wall_time_s'],
                    'Max_RSS_MB': res['max_rss_mb'], 'CPU_pct': res['cpu_pct'],
                    'IO_In_MB': res['io_in_mb'], 'IO_Out_MB': res['io_out_mb'],
                    'Source_File': os.path.relpath(logfile, sdir),
                    'Source_Format': res['source_format'],
                })
                rows.append(row)
                stats_fmt[res['source_format'] or 'missing'] += 1

    if not rows:
        print("⚠️  未解析到任何资源日志！")
        return 1

    df = pd.DataFrame(rows)
    out_file = os.path.join(args.out_dir, "resource_summary.tsv")
    df.to_csv(out_file, sep='\t', index=False)

    print(f"✅ resource_summary.tsv: {len(df)} 行, {df['Tool'].nunique()} 个工具")
    print(f"   格式统计: time-f={stats_fmt['time-f']}, time-v={stats_fmt['time-v']}, missing={stats_fmt['missing']}")
    for tool in sorted(df['Tool'].unique()):
        tool_df = df[df['Tool'] == tool]
        print(f"   {tool}: {len(tool_df)} 条, "
              f"avg_time={tool_df['Wall_Time_s'].mean():.0f}s, "
              f"avg_rss={tool_df['Max_RSS_MB'].mean():.0f}MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
