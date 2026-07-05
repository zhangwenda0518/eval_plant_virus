#!/usr/bin/env python3
"""
【宏基因组评估神器】 TSV 汇总报告转 CAMI 格式工具 (终极完美合并版)
特性:
1. 自动识别表头所有软件 (Pct_xxx)
2. 多进程极速批量处理
3. 严格遵循 OPAL/CAMI 格式规范 (精确生成 8 层级 TAXPATH)
4. 自动净化 SampleID，剥离无用后缀 (如 _abundance_PE, _PE) 与 Gold Standard 对齐
5. 【核心修复】完美兼容 NCBI 病毒界 (Viruses) 识别，并强制阻断伪 Root (TaxID 1) 污染。
"""

import os
import csv
import argparse
import pathlib
from concurrent.futures import ProcessPoolExecutor

# CAMI 官方的 8 大标准级别
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

# 映射表：将各软件的缩写映射为标准全拼
# 修复说明：
# 1. 恢复了 R 和 R1，因为很多软件（如 Kraken/Centrifuger）会将 Viruses (10239) 标记为 Root 的子节点 (R1)
# 2. 伪 Root (TaxID 1) 会在下方的代码逻辑中被直接拦截，确保不会被错误写成 superkingdom。
RANK_CODES = {
    'D': 'superkingdom', 'K': 'superkingdom', 'd': 'superkingdom', 'k': 'superkingdom',
    'R': 'superkingdom', 'R1': 'superkingdom',  
    'P': 'phylum', 'p': 'phylum',
    'C': 'class', 'c': 'class',
    'O': 'order', 'o': 'order',
    'F': 'family', 'f': 'family',
    'G': 'genus', 'g': 'genus',
    'S': 'species', 's': 'species',
    'S1': 'strain', 'S2': 'strain', 't': 'strain'
}

def process_single_file(input_file: str, output_dir: str):
    filename = os.path.basename(input_file)
    
    # 1. 剥离基本扩展名
    raw_id = filename.replace('_tree_report.tsv', '').replace('.tsv', '')
    # 2. 剥离多余的后缀，恢复纯净名称
    sample_id = raw_id.replace("_abundance_PE", "").replace("_PE", "")
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        if not lines: 
            return
            
        header = lines[0].strip('\n').split('\t')
        tools = [col.replace('Pct_', '') for col in header if col.startswith('Pct_')]
        if not tools: 
            return

        taxid_idx = header.index('TaxID')
        rank_idx = header.index('Rank')
        name_idx = header.index('Indented_Name')
        pct_indices = {tool: header.index(f'Pct_{tool}') for tool in tools}

        for tool in tools:
            output_filepath = os.path.join(output_dir, f"{sample_id}_{tool}.profile")
            with open(output_filepath, 'w', newline='', encoding='utf-8') as out_f:
                # 写入 CAMI 文件头 (此时写入的已经是净化后的 sample_id)
                out_f.write(f"@SampleID:\t{sample_id}\n")
                out_f.write(f"@Version:\t0.9.1\n")
                out_f.write(f"@Ranks:\t{'|'.join(RANKS)}\n")
                out_f.write(f"@ToolID:\t{tool}\n")
                
                writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')
                writer.writerow(["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"])
                
                # 状态机：用空字符串 "" 占位（绝不能用 "NA"）
                rank_state = {rank: {"taxid": "", "name": ""} for rank in RANKS}
                stack = []
                
                for line in lines[1:]:
                    if not line.strip(): continue
                    cols = line.strip('\n').split('\t')
                    
                    raw_name = cols[name_idx]
                    taxid = cols[taxid_idx]
                    raw_rank = cols[rank_idx]
                    
                    try:
                        pct = float(cols[pct_indices[tool]])
                    except ValueError: 
                        continue
                    if pct <= 0: 
                        continue 
                        
                    clean_name = raw_name.strip()
                    # 依据前导空格数判断进化树层级深度
                    indent_spaces = len(raw_name) - len(raw_name.lstrip(' '))
                    depth = max(0, indent_spaces // 2)
                    
                    # 弹出比当前深度更深的祖先节点
                    while len(stack) > depth:
                        popped_rank, _, _ = stack.pop()
                        mapped = RANK_CODES.get(popped_rank)
                        if mapped in rank_state:
                            rank_state[mapped] = {"taxid": "", "name": ""}
                            
                    stack.append((raw_rank, taxid, clean_name))
                    mapped_rank = RANK_CODES.get(raw_rank)
                    
                    if mapped_rank:
                        # 记录当前级别信息
                        rank_state[mapped_rank] = {"taxid": taxid, "name": clean_name}
                        # 清空当前级别以下的所有子级别
                        for lower in RANKS[RANKS.index(mapped_rank) + 1:]:
                            rank_state[lower] = {"taxid": "", "name": ""}
                            
                    # 【核心修复】：强制阻断 Root (1)，仅输出合法的生物分类 TaxID
                    if mapped_rank and taxid not in {"0", "1", "", "NA"}:
                        # 强制组装 8 个元素，严格遵守 OPAL 的 7 竖线标准
                        taxpath = "|".join([rank_state[r]["taxid"] for r in RANKS])
                        taxpathsn = "|".join([rank_state[r]["name"] for r in RANKS])
                        writer.writerow([taxid, mapped_rank, taxpath, taxpathsn, f"{pct:.6f}"])
    except Exception as e:
        print(f"[错误] 处理 {input_file} 时发生异常: {e}")

def main():
    parser = argparse.ArgumentParser(description="Convert Metagenomic TSV to CAMI profiling format.")
    parser.add_argument('-i', '--input', nargs='+', required=True, help="Input TSV files (supports wildcards)")
    parser.add_argument('-o', '--outdir', default='cami_profiles', help="Output directory")
    parser.add_argument('-t', '--threads', type=int, default=min(8, os.cpu_count() or 1), help="Number of threads")
    args = parser.parse_args()

    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    input_files = [os.path.abspath(f) for f in args.input if os.path.isfile(f)]
    
    print(f"🚀 开始转换并净化 {len(input_files)} 个文件 (线程数: {args.threads})...")
    
    if args.threads > 1 and len(input_files) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            for f in input_files: 
                executor.submit(process_single_file, f, args.outdir)
    else:
        for f in input_files: 
            process_single_file(f, args.outdir)
            
    print(f"✅ 完美！CAMI 评估文件已全部转换并完成 SampleID 与 Root 节点修正，保存在: {os.path.abspath(args.outdir)}")

if __name__ == "__main__":
    main()
