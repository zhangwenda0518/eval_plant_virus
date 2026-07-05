#!/usr/bin/env python3
"""
【宏基因组评估神器】 自动生成 CAMI 金标准 (Gold Standard) 智能通用版
特性：
1. 支持标准 5 列输入格式 (Sample_Name, Acc_num, TaxID, Reads, Abundance)
2. 【核心升级】如果缺少 Abundance 列，将自动根据本样本内的 Reads 数量计算相对丰度 (%)！
3. 自动填补合并单元格 (如 Excel 中复制出的空格)。
4. 强制毒株级别 (Strain) 降维补丁。
5. 【核心修复】强制补齐缺失的病毒界 (Viruses, 10239) 节点，彻底绝杀 OPAL 评估中的 nan 错误。
"""

import os
import csv
import argparse
import sys
import pandas as pd
from Bio import Entrez
from urllib.error import HTTPError, URLError

# CAMI 官方的 8 大标准级别
TARGET_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

def get_col(col_input):
    """辅助函数：如果是数字字符串则转为整型索引，否则保持列名"""
    if isinstance(col_input, str) and col_input.lstrip('-').isdigit():
        return int(col_input)
    return col_input

def generate_gold_standards(input_file, out_dir, email, c_samp, c_tax, c_reads, c_pct):
    if not os.path.exists(input_file):
        print(f"[错误] 找不到输入文件: {input_file}")
        sys.exit(1)
        
    print(f"📥 正在读取数据文件: {input_file} ...")
    try:
        # 使用 pandas 读取 TSV/CSV，不假设表头存在，允许通过数字索引
        sep = ',' if input_file.endswith('.csv') else '\t'
        # 如果列索引是纯数字，我们不把第一行当作表头，防止漏读数据
        header_row = None if str(c_samp).isdigit() else 0
        df = pd.read_csv(input_file, sep=sep, header=header_row, skipinitialspace=True)
        
        c_samp, c_tax, c_reads, c_pct = get_col(c_samp), get_col(c_tax), get_col(c_reads), get_col(c_pct)
        
        # ⚠️ 填补像 Excel 那样第一列留空的合并单元格
        df[c_samp] = df[c_samp].replace(r'^\s*$', pd.NA, regex=True).ffill()
        
        # 去除非数据行（如果第一行是表头文字，转换为数字时会变成 NaN，将其剔除）
        df = df[pd.to_numeric(df[c_tax], errors='coerce').notnull()].copy()
        
        # 格式化基础列
        df['Sample'] = df[c_samp].astype(str).str.strip().str.replace(" ", "")
        df['TaxID'] = pd.to_numeric(df[c_tax], errors='coerce').astype(int).astype(str)
        df['Reads'] = pd.to_numeric(df[c_reads], errors='coerce').fillna(0)
        
        # 🌟 【核心逻辑】判断是否需要根据 Reads 计算丰度
        need_calc = False
        if c_pct == -1: # 用户指定不需要/没有丰度列
            need_calc = True
        else:
            try:
                df['Abundance'] = pd.to_numeric(df[c_pct], errors='coerce')
                # 如果这一列全是 NaN，说明空有其表
                if df['Abundance'].isna().all():
                    need_calc = True
            except (KeyError, IndexError):
                need_calc = True
                
        if need_calc:
            print(f"💡 [智能提示] 未检测到有效的丰度列数据，正在根据 Reads 列自动计算每个样本内的相对丰度 (%)...")
            # 根据 Sample 分组，计算占比并乘以 100
            total_reads_per_sample = df.groupby('Sample')['Reads'].transform('sum')
            df['Abundance'] = (df['Reads'] / total_reads_per_sample) * 100
            # 填补可能出现的除以 0 的情况
            df['Abundance'] = df['Abundance'].fillna(0)
        else:
            print(f"💡 [提示] 成功读取丰度 (Abundance) 列数据。")

    except Exception as e:
        print(f"[错误] 解析表格失败: {e}\n请使用 -h 参数查看列索引的指定方法。")
        sys.exit(1)

    # 按 Sample 和 TaxID 汇总丰度 (合并同一个 Dataset 里相同的 TaxID)
    grouped = df.groupby(['Sample', 'TaxID'])['Abundance'].sum().reset_index()
    
    datasets = {}
    for _, row in grouped.iterrows():
        sample_name = str(row['Sample'])
        if sample_name not in datasets:
            datasets[sample_name] = {}
        datasets[sample_name][row['TaxID']] = row['Abundance']

    all_taxids = list(df['TaxID'].unique())
    print(f"📊 解析完毕: 发现 {len(datasets)} 个样本，共计 {len(all_taxids)} 个唯一 TaxID。")

    # 通过 NCBI 获取分类树
    Entrez.email = email
    print(f"📡 正在通过 NCBI Entrez 接口查询分类树 (这可能需要几秒钟)...")
    try:
        handle = Entrez.efetch(db="taxonomy", id=",".join(all_taxids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except (HTTPError, URLError) as e:
        print(f"[错误] NCBI 网络请求失败: {e}。请检查网络，或稍后再试。")
        sys.exit(1)
    except Exception as e:
        print(f"[错误] NCBI 接口异常: {e}")
        sys.exit(1)

    tax_lineages = {}
    for rec in records:
        tid = rec["TaxId"]
        lineage = rec["LineageEx"]
        lineage.append({'TaxId': tid, 'Rank': rec['Rank'], 'ScientificName': rec['ScientificName']})
        tax_lineages[tid] = lineage

    print("✅ NCBI 分类树获取成功！")

    # 生成 CAMI 格式文件
    os.makedirs(out_dir, exist_ok=True)
    generated_count = 0

    for ds_name, taxid_pcts in datasets.items():
        profile_path = os.path.join(out_dir, f"{ds_name}_GoldStandard.profile")
        rank_abundances = {}
        
        for tid, pct in taxid_pcts.items():
            if tid not in tax_lineages: 
                print(f"  [警告] NCBI 中找不到 TaxID {tid}，忽略该项。")
                continue
                
            lineage = tax_lineages[tid]
            path_ids = [""] * 8
            path_names = [""] * 8
            
            for node in lineage:
                r = node.get('Rank', '')
                if r in TARGET_RANKS:
                    idx = TARGET_RANKS.index(r)
                    path_ids[idx] = node['TaxId']
                    path_names[idx] = node['ScientificName']
                    
            # 【毒株级特供补丁】强制将底层 TaxID 放入 Strain 级别
            if tid not in path_ids:
                path_ids[7] = tid
                for node in lineage:
                    if node['TaxId'] == tid:
                        path_names[7] = node['ScientificName']
                        break
                        
            # 🌟🌟🌟【核心修复：病毒界强制填补补丁】🌟🌟🌟
            # 如果 NCBI 分类中 superkingdom 为空，则强制指定为 Viruses (10239)
            if path_ids[0] == "":
                path_ids[0] = "10239"
                path_names[0] = "Viruses"
                        
            for i in range(8):
                if path_ids[i] != "":
                    # 截取到当前级别，后面的层级强制用 "" 补齐
                    current_path_ids = path_ids[:i+1] + [""] * (7 - i)
                    current_path_names = path_names[:i+1] + [""] * (7 - i)
                    path_key = "|".join(current_path_ids)
                    
                    if path_key not in rank_abundances:
                        rank_abundances[path_key] = {
                            "taxid": path_ids[i], "rank": TARGET_RANKS[i],
                            "taxpath": path_key, "taxpathsn": "|".join(current_path_names), "pct": 0.0
                        }
                    rank_abundances[path_key]["pct"] += pct

        # 写入文件
        with open(profile_path, 'w', encoding='utf-8') as out_f:
            out_f.write(f"@SampleID:\t{ds_name}\n@Version:\t0.9.1\n")
            out_f.write(f"@Ranks:\t{'|'.join(TARGET_RANKS)}\n@ToolID:\tGoldStandard\n")
            writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')
            writer.writerow(["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"])
            for node in rank_abundances.values():
                writer.writerow([node["taxid"], node["rank"], node["taxpath"], node["taxpathsn"], f"{node['pct']:.6f}"])
                
        generated_count += 1
                
    print(f"\n🎉 完美！成功生成 {generated_count} 个金标准文件，保存在: {os.path.abspath(out_dir)}")

def main():
    example_text = '''
【输入文件格式说明】
输入文件通常为一个 5 列的 TSV 或 CSV 文件：
Sample_Name (样本名) | Acc_num (登记号) | TaxID | Reads (读长数) | Abundance (丰度)

如果你的表格中没有 Abundance 这一列，或者这一列为空，程序将智能地
根据本样本内的 Reads 数量，自动为你计算出百分比相对丰度 (Abundance %)！

[例子 1：标准的 5 列文件]
假设列索引分别为 0, 1, 2, 3, 4 (索引从 0 开始)。
命令: python make_gold.py -i truth.tsv -s 0 -t 2 -r 3 -p 4

[例子 2：只有 4 列的文件，缺少丰度列]
你可以将丰度列索引设为 -1，程序将直接利用第 3 列 (Reads) 自动计算丰度。
命令: python make_gold.py -i truth.tsv -s 0 -t 2 -r 3 -p -1
    '''

    parser = argparse.ArgumentParser(
        description="向 NCBI 查询分类树，生成 CAMI 金标准文件。支持 Reads 丰度自动计算。",
        epilog=example_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True, help="包含 Ground Truth 的输入表格文件 (TSV/CSV)。")
    parser.add_argument('-o', '--outdir', default="gold_standards", help="输出文件夹的名称 (默认: gold_standards)")
    parser.add_argument('-e', '--email', default="benchmark@example.com", help="NCBI API 要求的联络邮箱")
    
    # 列索引参数
    parser.add_argument('-s', '--sample-col', default=0, help="样本名称 (Sample_Name) 所在的列索引 (默认: 0)")
    parser.add_argument('-t', '--taxid-col', default=2, help="TaxID 所在的列索引 (默认: 2)")
    parser.add_argument('-r', '--reads-col', default=3, help="Reads 数量所在的列索引 (默认: 3)")
    parser.add_argument('-p', '--pct-col', default=4, help="丰度 (Abundance) 所在的列索引。如果缺失，设为 -1 自动计算 (默认: 4)")

    args = parser.parse_args()
    generate_gold_standards(args.input, args.outdir, args.email, args.sample_col, args.taxid_col, args.reads_col, args.pct_col)

if __name__ == "__main__":
    main()
