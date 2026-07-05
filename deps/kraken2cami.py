#!/usr/bin/env python3
"""
【宏基因组评估神器】 Kraken2 原生报告转 CAMI 格式工具 (终极批量版)
用途: 将 Kraken2 或 Bracken 生成的标准 6 列 report.tsv 批量转换为严格的 CAMI 格式。
特性:
1. 完美兼容 OPAL (消除 Index/Key Error，强制 8 级空字符串对齐)。
2. 支持毒株级别 (Strain/S1/S2)。
3. 多进程并发处理，自动提取 SampleID。
"""

import os
import csv
import argparse
import pathlib
from concurrent.futures import ProcessPoolExecutor

# ==========================================
# 1. 核心配置与映射字典
# ==========================================
# CAMI 官方的 8 大标准级别
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

# 映射表：将 Kraken2 简写映射为标准全拼
RANK_CODES = {
    'D': 'superkingdom', 'R1': 'superkingdom', 'R': 'superkingdom', 'K': 'superkingdom',
    'P': 'phylum', 'C': 'class', 'O': 'order', 'F': 'family',
    'G': 'genus', 'S': 'species', 'S1': 'strain', 'S2': 'strain'
}

# ==========================================
# 2. 核心转换逻辑
# ==========================================
def process_single_kreport(input_file: str, output_dir: str, tool_name: str):
    """解析单个 Kraken2 报告并生成 CAMI profile"""
    # 自动从文件名提取 SampleID (去除常见的后缀)
    filename = os.path.basename(input_file)
    sample_id = filename.replace('.kreport', '').replace('.report', '').replace('.tsv', '').replace('.txt', '')
    output_filepath = os.path.join(output_dir, f"{sample_id}_{tool_name}.profile")
    
    try:
        with open(input_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            
        if not lines:
            print(f"[警告] 文件为空: {input_file}")
            return

        with open(output_filepath, "w", newline="", encoding="utf-8") as out_f:
            # 写入标准的 CAMI Header
            out_f.write(f"@SampleID:\t{sample_id}\n")
            out_f.write(f"@Version:\t0.9.1\n")
            out_f.write(f"@Ranks:\t{'|'.join(RANKS)}\n")
            out_f.write(f"@ToolID:\t{tool_name}\n")
            
            writer = csv.writer(out_f, delimiter="\t", lineterminator="\n")
            writer.writerow(["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"])
            
            # 状态机：默认空字符串 ""
            rank_state = {rank: {"taxid": "", "name": ""} for rank in RANKS}
            stack = []
            
            for line in lines:
                if not line.strip(): continue
                parts = line.rstrip("\n").split("\t")
                
                # 严格假设这里是 Kraken2 的 6列原生格式!
                if len(parts) < 6: continue
                
                try:
                    perc = float(parts[0]) # 第 0 列是丰度百分比
                except ValueError:
                    continue
                if perc <= 0.0:
                    continue
                    
                raw_rank_code = parts[3].strip() # 第 3 列是 Rank 层级 (如 D, P, C, S, S1)
                if not raw_rank_code: continue
                
                # 支持 S1 这样的双字符 Rank，回退机制为取首字母
                rank_code = raw_rank_code.upper() 
                if rank_code not in RANK_CODES and len(rank_code) > 0:
                    rank_code = rank_code[0] 
                    
                taxid = parts[4].strip()     # 第 4 列是 TaxID
                name_field = parts[5]        # 第 5 列是 Indented Name
                
                # 计算缩进以判断深度 (Kraken2 默认每层缩进 2 个空格)
                leading_spaces = len(name_field) - len(name_field.lstrip())
                depth = max(0, leading_spaces // 2)
                clean_name = name_field.strip() or taxid or "unknown"

                # 弹出比当前深度更深的节点
                while len(stack) > depth:
                    popped_code, _, _ = stack.pop()
                    mapped = RANK_CODES.get(popped_code)
                    if mapped in rank_state:
                        rank_state[mapped] = {"taxid": "", "name": ""}

                stack.append((rank_code, taxid, clean_name))
                mapped_rank = RANK_CODES.get(rank_code)

                if mapped_rank:
                    # 记录当前级别信息
                    rank_state[mapped_rank] = {"taxid": taxid, "name": clean_name}
                    # 清空当前级别以下的所有子级别
                    for lower in RANKS[RANKS.index(mapped_rank) + 1:]:
                        rank_state[lower] = {"taxid": "", "name": ""}

                # 如果它是标准 8 个级别之一，输出该行
                if mapped_rank and taxid not in {"0", "", "NA"}:
                    # 强制组装 8 个元素，严格遵守 OPAL 的 7 竖线标准
                    taxpath = "|".join([rank_state[r]["taxid"] for r in RANKS])
                    taxpathsn = "|".join([rank_state[r]["name"] for r in RANKS])
                    writer.writerow([taxid, mapped_rank, taxpath, taxpathsn, f"{perc:.6f}"])
                    
    except Exception as e:
        print(f"[错误] 处理 {input_file} 时发生异常: {e}")

# ==========================================
# 3. 命令行接口与多进程调度
# ==========================================
def main():
    example_text = '''
【使用示例】
1. 处理单个文件:
   python kraken2cami.py -i sample1.kreport -o output_dir --tool kraken2

2. 批量处理当前目录下所有的 Kraken2 报告 (开启 8 线程极速转换):
   python kraken2cami.py -i *.kreport -o output_dir -t 8 --tool kraken2
    '''
    
    parser = argparse.ArgumentParser(
        description="将 Kraken2/Bracken 的标准 6 列 report 报告批量转换为 CAMI 格式。",
        epilog=example_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', nargs='+', required=True, help="输入的 Kraken2 报告文件 (支持通配符 *.kreport)")
    parser.add_argument('-o', '--outdir', default='cami_profiles', help="输出目录 (默认: cami_profiles)")
    parser.add_argument('--tool', default='kraken2', help="工具名称，将写入 CAMI 文件头中 (默认: kraken2)")
    parser.add_argument('-t', '--threads', type=int, default=min(8, os.cpu_count() or 1), help="多进程线程数 (默认: 自动检测)")

    args = parser.parse_args()

    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    input_files = [os.path.abspath(f) for f in args.input if os.path.isfile(f)]
    
    if not input_files:
        print("[错误] 未找到任何有效的输入文件，请检查路径。")
        return

    print(f"🚀 开始转换 {len(input_files)} 个 Kraken2 报告 (线程数: {args.threads}, 标记为: {args.tool})...")
    
    if args.threads > 1 and len(input_files) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            for f in input_files: 
                executor.submit(process_single_kreport, f, args.outdir, args.tool)
    else:
        for f in input_files: 
            process_single_kreport(f, args.outdir, args.tool)
            
    print(f"✅ 转换完毕！完美的 CAMI 文件已保存在: {os.path.abspath(args.outdir)}")

if __name__ == "__main__":
    main()
