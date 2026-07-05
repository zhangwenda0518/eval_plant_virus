#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Virus Assembly Benchmark — 公共工具模块 (已针对深度梯度及节段病毒优化)
================================================================================
供 benchmark_mq.py / benchmark_chimeric.py / benchmark_resource.py / benchmark_summarize.py 共享使用。
"""

import os, re, subprocess, glob
from pathlib import Path

# ====================================================================
# 全局配置
# ====================================================================

# 10 种工具的配色 (colorblind-friendly)
MY_PALETTE = {
    "Megahit":              "#4C72B0",   # 蓝
    "RNAViralSPAdes":       "#55A868",   # 绿
    "Penguin":              "#DD8452",   # 橙
    "MetaViralSPAdes":      "#A0522D",   # 褐
    "Trinity":              "#6BAED6",   # 浅蓝
    "RNABloom":             "#E69F00",   # 金
    "MH_Merge":             "#8172B3",   # 紫
    "MH_Split_Merge":       "#F5A623",   # 深金
    "ALL_Merge":            "#DA8BC3",   # 粉
    "RefineC_Merge":        "#C44E52",   # 红
}

MUT_ORDER  = ["0.0", "5.0", "10.0", "15.0"]
BG_ORDER   = ["0.1", "0.5", "1", "2", "5", "10"]

# 🌟 更新：包含 5M 和 Master（以 10000000 代表）
DEPTH_ORDER = [10, 20, 50, 100, 200, 500, 1000, 2000]
# MetaQUAST 质量评估关键指标 (QUAST 标准输出列名)
METAQUAST_METRICS = [
    "Genome fraction (%)",
    "Largest alignment",
    "NGA50",
    "# mismatches per 100 kbp",
    "# indels per 100 kbp",
    "# misassemblies",
    "Duplication ratio",
    "# contigs",
    "Reference length",  # 用于计算相对 NGA50 / Largest alignment
]

# 评分权重 (Genome fraction + Largest alignment/RefLen 各 30%)
SCORE_WEIGHTS = {
    "Genome fraction (%)":         0.30,
    "Largest_alignment_relative":  0.30,
    "NGA50_relative":              0.10,
    "# mismatches per 100 kbp":    0.10,
    "# misassemblies":             0.10,
    "Duplication ratio":           0.10,
}

# 越大越好的指标
BIG_BETTER = {"Genome fraction (%)", "Largest alignment", "NGA50",
              "Largest_alignment_relative", "NGA50_relative"}

# 文件命名模式
FILE_PATTERNS = {
    "Megahit":        "{base}_megahit.contig.fasta",
    "RNAViralSPAdes": "{base}_rnaviralspades.contig.fasta",
    "Penguin":        "{base}_penguin.contig.fasta",
    "MetaViralSPAdes":"{base}_metaviralspades.contig.fasta",
    "Trinity":        "{base}_trinity.contig.fasta",
    "RNABloom":       "{base}_rnabloom.contig.fasta",
    "MH_Merge":       "{base}_MH_merge.merged.fasta",
    "MH_Split_Merge": "{base}_MH_split_merge.merged.fasta",
    "ALL_Merge":      "{base}_ALL_merge.merged.fasta",
    "RefineC_Merge":  "{base}_all_tools_refineC_merge.merged.fasta",
}

BASE_TOOLS = ["Megahit", "RNAViralSPAdes", "Penguin", "MetaViralSPAdes", "Trinity", "RNABloom"]

MERGE_TOOL_SUBDIR = {
    "MH_Merge":       "MH_merge",
    "MH_Split_Merge": "MH_split_merge",
    "ALL_Merge":      "ALL_merge",
}

REFINEC_SUBDIR_SUFFIX = "_all_tools_refineC_merge"


# ====================================================================
# 工具函数
# ====================================================================

def get_existing_file(filepath):
    """匹配 .fasta 或 .fasta.gz"""
    for p in [filepath, filepath + ".gz"]:
        if os.path.exists(p):
            return p
    return None

def parse_dir_name(dirname):
    """提取基础变量（兼容旧版、sub 以及 depth 格式）"""
    result = {}
    m = re.search(r'(\d+\.?\d*)pct', dirname)
    if m: result['mutation_rate'] = m.group(1)
    
    # 🌟 全面兼容 sub 和 depth 格式
    base_lower = dirname.lower()
    if "master" in base_lower:
        result['depth'] = 10000000
    else:
        # 支持 _sub_100, _depth1000, sub_1M 等等
        m = re.search(r'(?:sub|depth)_?(\d+\.?\d*[kKmM]?)', dirname, re.IGNORECASE)
        if m:
            val = m.group(1).upper()
            if val.endswith('K'):
                result['depth'] = int(float(val[:-1]) * 1000)
            elif val.endswith('M'):
                result['depth'] = int(float(val[:-1]) * 1000000)
            else:
                result['depth'] = int(val)
        else:
            result['depth'] = 0

    m = re.search(r'_rep(\d+)', dirname)
    if m: result['rep'] = f"rep{m.group(1)}"
    m = re.search(r'group(\d+)', dirname)
    if m: result['group'] = int(m.group(1))
    # 背景比例: bg_r0.1, bg_r1, bg_r10 等
    m = re.search(r'bg_r(\d+\.?\d*)', dirname, re.IGNORECASE)
    if m: result['background_ratio'] = m.group(1)
    return result


def parse_resource_log(filepath):
    """解析资源消耗日志 (双格式兼容)"""
    result = {
        'wall_time_s': None, 'max_rss_mb': None, 'cpu_pct': None,
        'io_in_mb': None, 'io_out_mb': None, 'source_format': None
    }
    if filepath is None or not os.path.exists(filepath):
        return result
    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except Exception:
        return result

    if 'Elapsed (wall clock) time' in content:
        result['source_format'] = 'time-v'
        m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', content)
        if m:
            parts = m.group(1).strip().split(':')
            if len(parts) == 3:
                result['wall_time_s'] = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
            elif len(parts) == 2:
                result['wall_time_s'] = float(parts[0]) * 60 + float(parts[1])
            else:
                result['wall_time_s'] = float(parts[0])
        m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', content)
        if m: result['max_rss_mb'] = float(m.group(1)) / 1024.0
        m = re.search(r'Percent of CPU this job got:\s*(\d+)%', content)
        if m: result['cpu_pct'] = float(m.group(1))
        m = re.search(r'File system inputs:\s*(\d+)', content)
        if m: result['io_in_mb'] = float(m.group(1)) * 512 / (1024**2)
        m = re.search(r'File system outputs:\s*(\d+)', content)
        if m: result['io_out_mb'] = float(m.group(1)) * 512 / (1024**2)

    elif 'Time:' in content and 'Memory:' in content:
        result['source_format'] = 'time-f'
        m = re.search(r'Time:\s*([\d.]+)\s*seconds', content)
        if m: result['wall_time_s'] = float(m.group(1))
        m = re.search(r'Memory:\s*(\d+)\s*KB', content)
        if m: result['max_rss_mb'] = float(m.group(1)) / 1024.0
        m = re.search(r'CPU:\s*(\d+)%', content)
        if m: result['cpu_pct'] = float(m.group(1))
    return result


def discover_samples(input_dir):
    """
    扫描 input_dir 下的样本目录。
    🌟 已经针对深度梯度模式及大小写格式完美适配。
    """
    samples = []
    # 模式 A: 旧版 (带突变率和分组)
    pattern_a = os.path.join(input_dir, "group*_mut*", "Master_*pct_sub_*_rep*")
    for sdir in sorted(glob.glob(pattern_a)):
        if not os.path.isdir(sdir): continue
        base_name = os.path.basename(os.path.normpath(sdir))
        info = parse_dir_name(base_name)
        parent = os.path.basename(os.path.dirname(os.path.normpath(sdir)))
        pinfo = parse_dir_name(parent)
        info.setdefault('group', pinfo.get('group', 0))
        info.setdefault('mutation_rate', pinfo.get('mutation_rate', '0.0'))
        
        # 为了兼容不同的调用脚本，双写大小写键
        info['Depth'] = info.get('depth', 0)
        info['Mutation_Rate'] = info['mutation_rate']
        info['Group'] = info['group']
        info['Rep'] = info.get('rep', 'rep1')
        
        info['sample_dir'] = sdir
        info['base_name'] = base_name
        samples.append(info)

    if samples:
        return samples

    # 🌟 模式 B: eval*_sub_* 和 eval*_master  (扁平化深度梯度专用)
    pattern_b = os.path.join(input_dir, "eval*")
    for sdir in sorted(glob.glob(pattern_b)):
        if not os.path.isdir(sdir):
            continue
        base_name = os.path.basename(os.path.normpath(sdir))
        info = {}
        
        # --- 🌟 智能解析深度 (适配 eval3_depthXXX_PE 格式) ---
        base_lower = base_name.lower()
        if "master" in base_lower:
            depth_val = 10000000
        else:
            # 🌟 注意这里：把 'sub_' 改成了 'depth'
            m = re.search(r'depth(\d+\.?\d*[kKmM]?)', base_name, re.IGNORECASE)
            if m:
                val = m.group(1).upper()
                if val.endswith('K'):
                    depth_val = int(float(val[:-1]) * 1000)
                elif val.endswith('M'):
                    depth_val = int(float(val[:-1]) * 1000000)
                else:
                    depth_val = int(val)
            else:
                depth_val = 0
                
        # 兼容大小写写入
        bg_ratio = parse_dir_name(base_name).get('background_ratio', None)
        info['depth'] = info['Depth'] = depth_val
        info['mutation_rate'] = info['Mutation_Rate'] = '0.0'
        info['background_ratio'] = info['Background_Ratio'] = bg_ratio
        m = re.search(r'_rep(\d+)', base_name)
        info['rep'] = info['Rep'] = f"rep{m.group(1)}" if m else 'rep1'
        info['group'] = info['Group'] = 0

        info['sample_dir'] = sdir
        info['base_name'] = base_name
        samples.append(info)

    if samples:
        return samples

    # 🌟 模式 C: 单病毒 LoD_mix 输出 {accession}_depth{N}_rep{R}/
    pattern_c = os.path.join(input_dir, "*_depth*_rep*")
    for sdir in sorted(glob.glob(pattern_c)):
        if not os.path.isdir(sdir):
            continue
        base_name = os.path.basename(os.path.normpath(sdir))
        info = {}
        m = re.search(r'depth(\d+)', base_name)
        depth_val = int(m.group(1)) if m else 0
        info['depth'] = info['Depth'] = depth_val
        info['mutation_rate'] = info['Mutation_Rate'] = '0.0'
        m = re.search(r'_rep(\d+)', base_name)
        info['rep'] = info['Rep'] = f"rep{m.group(1)}" if m else 'rep1'
        info['group'] = info['Group'] = 0
        info['sample_dir'] = sdir
        info['base_name'] = base_name
        samples.append(info)

    return samples


def get_assembly_files(sample_dir, base_name, mode='7'):
    if mode in ('7', '10'):
        active_tools = list(FILE_PATTERNS.keys())
    elif mode == '6':
        active_tools = ["Megahit", "RNAViralSPAdes", "Penguin", "MetaViralSPAdes", "Trinity", "RNABloom"]
    else:
        active_tools = ["Megahit", "RNAViralSPAdes", "Penguin", "RefineC_Merge"]
    results = []
    for tool in active_tools:
        expected = os.path.join(sample_dir, FILE_PATTERNS[tool].format(base=base_name))
        path = get_existing_file(expected)
        if path:
            results.append((tool, path))
    return results
