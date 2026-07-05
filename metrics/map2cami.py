#!/usr/bin/env python3
"""
【宏基因组评估】 V40 all_viruses.best.summary.tsv 转 CAMI 工具 (自带 SampleID 智能净化版)
特性:
1. 专为 V40 病毒定量全景文件设计，完美适配多维度指标。
2. 【核心提升】集成 SampleID 自动净化功能，默认剔除 _abundance_PE、_PE 等后缀，无损对齐 Gold Standard！
3. 支持高阶正则表达式清洗，从容应对各类测序批次尾缀。
"""

import os
import csv
import argparse
import sys
import re
import pandas as pd
from Bio import Entrez

TARGET_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]

def clean_sample_name(raw_name: str, regex_pattern: str = None) -> str:
    """样本名净化核心函数"""
    sample_id = str(raw_name).strip()
    
    # 1. 基础硬编码清洗（集成 fix_sampleid.py 的核心逻辑）
    sample_id = sample_id.replace("_abundance_PE", "").replace("_PE", "").replace("_SE", "")
    
    # 2. 如果用户传入了自定义的正则规则，执行高级清洗
    if regex_pattern:
        try:
            sample_id = re.sub(regex_pattern, "", sample_id, flags=re.IGNORECASE)
        except Exception as e:
            print(f"[警告] 正则表达式替换失败: {e}")
            
    # 去除首尾可能残余的下划线、连字符或空格
    return sample_id.strip("_- ")

def convert_v40_to_cami(input_file, out_dir, tool_name, email, c_samp, c_tax, c_value, gt_file, gt_samp_col, gt_reads_col, clean_regex):
    if not os.path.exists(input_file):
        print(f"[错误] 找不到输入文件: {input_file}")
        sys.exit(1)
        
    print(f"📥 正在读取并解析 V40 病毒定量全景文件: {input_file} ...")
    try:
        # 显式规定使用 \t 分隔符
        df = pd.read_csv(input_file, sep='\t', skipinitialspace=True)
        
        # 检查关键列
        for col in [c_samp, c_tax, c_value]:
            if col not in df.columns:
                print(f"[错误] 文件中未找到指定的列: '{col}'。")
                print(f"当前文件包含的列有: {', '.join(df.columns.tolist()[:8])} ...")
                sys.exit(1)
                
        print(f"   => 成功映射列: 样本->[{c_samp}], TaxID->[{c_tax}], 定量指标->[{c_value}]")
        
        # 提取核心列
        df_core = df[[c_samp, c_tax, c_value]].copy()
        df_core.columns = ['Raw_Sample', 'TaxID', 'Value']
        
        # 🌟 【核心提升】在这里直接清洗 SampleID
        df_core['Sample'] = df_core['Raw_Sample'].apply(lambda x: clean_sample_name(x, clean_regex))
        
        # 打印部分清洗日志，让用户心里有数
        sample_mapping_example = df_core[['Raw_Sample', 'Sample']].drop_duplicates().head(3)
        print("💡 [SampleID 净化预览]：")
        for _, row in sample_mapping_example.iterrows():
            print(f"   => 原名称: {row['Raw_Sample']} ----> 净化后: {row['Sample']}")
            
        df_core = df_core[pd.to_numeric(df_core['TaxID'], errors='coerce').notnull()]
        df_core['TaxID'] = df_core['TaxID'].astype(int).astype(str)
        df_core['Value'] = pd.to_numeric(df_core['Value'], errors='coerce').fillna(0)
        
        # 🌟 核心丰度计算逻辑 🌟
        if gt_file:
            print("💡 [金标准对齐模式] 正在结合外部元数据重新校准真实丰度分母...")
            sep = ',' if gt_file.endswith('.csv') else '\t'
            df_gt = pd.read_csv(gt_file, sep=sep, skipinitialspace=True)
            g_s = df_gt.columns[int(gt_samp_col)] if str(gt_samp_col).isdigit() else gt_samp_col
            g_r = df_gt.columns[int(gt_reads_col)] if str(gt_reads_col).isdigit() else gt_reads_col
            
            df_gt['Sample_Key'] = df_gt[g_s].astype(str).str.strip().str.replace(" ", "")
            df_gt['True_Reads'] = pd.to_numeric(df_gt[g_r], errors='coerce').fillna(0)
            gt_dict = df_gt.groupby('Sample_Key')['True_Reads'].sum().to_dict()
            
            df_core['TrueTotal'] = df_core['Sample'].map(gt_dict)
            df_core['Abundance'] = (df_core['Value'] / df_core['TrueTotal']) * 100
            df_core['Abundance'] = df_core['Abundance'].fillna((df_core['Value'] / df_core.groupby('Sample')['Value'].transform('sum')) * 100)
        else:
            if 'abund' in c_value.lower():
                print("💡 [相对丰度直读] 直接提取 V40 自带的相对百分比 Asm_Rel_Abund(%)。")
                df_core['Abundance'] = df_core['Value']
            elif 'rpm' in c_value.lower():
                print("💡 [RPM全局标准化] 正在自动转换为以全局测序量为分母的绝对百分比(%%)。")
                df_core['Abundance'] = df_core['Value'] / 10000.0
            else:
                print("💡 [常规计数模式] 以当前样本内总 Mapped Reads 作为分母计算百分比。")
                df_core['Abundance'] = (df_core['Value'] / df_core.groupby('Sample')['Value'].transform('sum')) * 100
                
        df_core['Abundance'] = df_core['Abundance'].fillna(0)

    except Exception as e:
        print(f"[错误] 数据流解析失败: {e}")
        sys.exit(1)

    # 聚合去重
    grouped = df_core.groupby(['Sample', 'TaxID'])['Abundance'].sum().reset_index()
    datasets = {}
    for _, row in grouped.iterrows():
        s = str(row['Sample'])
        if s not in datasets: datasets[s] = {}
        datasets[s][row['TaxID']] = row['Abundance']

    all_taxids = list(df_core['TaxID'].unique())
    print(f"📊 过滤与汇总完毕：共 {len(datasets)} 个样本，包含 {len(all_taxids)} 个唯一 TaxID。")

    # 联网请求 NCBI Taxonomy
    Entrez.email = email
    print(f"📡 正在请求 NCBI 分类树回溯上级单元...")
    try:
        records = []
        batch_size = 300
        for i in range(0, len(all_taxids), batch_size):
            batch_ids = all_taxids[i:i+batch_size]
            handle = Entrez.efetch(db="taxonomy", id=",".join(batch_ids), retmode="xml")
            records.extend(Entrez.read(handle))
            handle.close()
    except Exception as e:
        print(f"[错误] NCBI 接口连接失败: {e}")
        sys.exit(1)

    tax_lineages = {rec["TaxId"]: rec["LineageEx"] + [{'TaxId': rec["TaxId"], 'Rank': rec['Rank'], 'ScientificName': rec['ScientificName']}] for rec in records}
    
    os.makedirs(out_dir, exist_ok=True)
    generated_count = 0
    
    for ds_name, taxid_pcts in datasets.items():
        profile_path = os.path.join(out_dir, f"{ds_name}_{tool_name}.profile")
        rank_abundances = {}
        
        for tid, pct in taxid_pcts.items():
            if tid not in tax_lineages: continue
            lineage = tax_lineages[tid]
            
            path_ids, path_names = [""] * 8, [""] * 8
            for node in lineage:
                r = node.get('Rank', '')
                if r in TARGET_RANKS:
                    idx = TARGET_RANKS.index(r)
                    path_ids[idx], path_names[idx] = node['TaxId'], node['ScientificName']
                    
            if tid not in path_ids:
                path_ids[7] = tid
                for node in lineage:
                    if node['TaxId'] == tid:
                        path_names[7] = node['ScientificName']
                        break
                        
            for i in range(8):
                if path_ids[i] != "":
                    current_path_ids = path_ids[:i+1] + [""] * (7 - i)
                    current_path_names = path_names[:i+1] + [""] * (7 - i)
                    path_key = "|".join(current_path_ids)
                    if path_key not in rank_abundances:
                        rank_abundances[path_key] = {
                            "taxid": path_ids[i], "rank": TARGET_RANKS[i],
                            "taxpath": path_key, "taxpathsn": "|".join(current_path_names), "pct": 0.0
                        }
                    rank_abundances[path_key]["pct"] += pct

        with open(profile_path, 'w', encoding='utf-8') as out_f:
            out_f.write(f"@SampleID:\t{ds_name}\n@Version:\t0.9.1\n")
            out_f.write(f"@Ranks:\t{'|'.join(TARGET_RANKS)}\n@ToolID:\t{tool_name}\n")
            writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')
            writer.writerow(["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"])
            for node in rank_abundances.values():
                writer.writerow([node["taxid"], node["rank"], node["taxpath"], node["taxpathsn"], f"{node['pct']:.6f}"])
        generated_count += 1
                
    print(f"✅ 生成完毕！{generated_count} 个净化后的 CAMI 预测文件已成功保存在: {os.path.abspath(out_dir)}")

def main():
    parser = argparse.ArgumentParser(description="V40 病毒定量结果转 CAMI 标准格式工具 (集成 SampleID 智能清洗)")
    parser.add_argument('-i', '--input', required=True, help="V40 输出的 all_viruses.best.summary.tsv 文件")
    parser.add_argument('-o', '--outdir', default="cami_profiles", help="输出目录")
    parser.add_argument('--tool', required=True, help="软件/流程名称")
    parser.add_argument('-e', '--email', default="your_email@example.com", help="NCBI Entrez 邮箱")
    
    # 针对 V40 优化的默认列名
    parser.add_argument('-s', '--sample-col', default="Sample", help="样本名列")
    parser.add_argument('-t', '--taxid-col', default="taxid", help="TaxID列 (V40中小写的 taxid)")
    parser.add_argument('-v', '--value-col', default="Asm_Rel_Abund(%)", help="定量列，可选 Asm_Rel_Abund(%%) 或 Asm_RPM 或 Asm_EM_Reads")
    
    # 【新增】SampleID 高级正则清理参数
    parser.add_argument('--clean-regex', default=None, help="高级自定义正则清洗规则（例如：输入 '(_lane\d+|_filtered)' 可以剔除附加的通道或过滤后缀）")
    
    # 外部金标准对齐（可选）
    parser.add_argument('--gt', help="外部金标准 Ground Truth 文件")
    parser.add_argument('--gt-samp', default="Sample", help="金标准样本列名")
    parser.add_argument('--gt-reads', default="Reads", help="金标准真实 Reads 数列")

    args = parser.parse_args()
    convert_v40_to_cami(args.input, args.outdir, args.tool, args.email, 
                        args.sample_col, args.taxid_col, args.value_col, 
                        args.gt, args.gt_samp, args.gt_reads, args.clean_regex)

if __name__ == "__main__":
    main()
