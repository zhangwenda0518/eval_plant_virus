#!/usr/bin/env python3
"""
【数据预处理】将丰度表与参考信息表合并，生成符合 make_cami.gold.py 要求的 5 列标准输入文件。
"""

import argparse
import sys
import pandas as pd

def merge_cami_data(gold_file, ref_file, out_file, to_pct):
    print(f"📥 正在读取丰度表: {gold_file}")
    # 使用 delim_whitespace=True 处理由多个空格或制表符分隔的文件
    df_gold = pd.read_csv(gold_file, delim_whitespace=True)

    print(f"📥 正在读取参考信息表: {ref_file}")
    # 读取包含 Taxid 的参考文件
    df_ref = pd.read_csv(ref_file, sep='\t', usecols=['Accession', 'Taxid'])

    print("🔄 正在合并数据并提取 Taxid...")
    df_merged = pd.merge(df_gold, df_ref, on='Accession', how='left')

    # 提取需要的 5 列
    try:
        df_final = df_merged[['Sample', 'Accession', 'Taxid', 'Expected_Reads', 'True_Abundance']].copy()
    except KeyError as e:
        print(f"❌ [错误] 缺失必要的列，请检查表头: {e}")
        sys.exit(1)

    # 检查是否有未匹配到的 Taxid
    missing_taxid = df_final['Taxid'].isna().sum()
    if missing_taxid > 0:
        print(f"⚠️ [警告] 发现 {missing_taxid} 条记录在 {ref_file} 中找不到对应的 Taxid！它们将被保留但 Taxid 为空。")

    # 处理丰度百分比转换
    if to_pct:
        print("🧮 正在将 True_Abundance 转换为百分比格式 (x100)...")
        df_final['True_Abundance'] = df_final['True_Abundance'] * 100

    print(f"💾 正在保存结果至: {out_file}")
    df_final.to_csv(out_file, sep='\t', index=False)
    print("✅ 处理完成！")

def main():
    parser = argparse.ArgumentParser(
        description="合并丰度表与参考信息表，生成 CAMI 金标准制作工具的输入文件。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # 必需参数
    parser.add_argument('-g', '--gold', required=True, help="[必需] 包含丰度和 Reads 信息的原始输入表 (例如: eval_known_virus_gold.tsv)")
    parser.add_argument('-r', '--ref', required=True, help="[必需] 包含 Accession 和 Taxid 对应关系的参考表 (例如: final.cluster.ref_info.tsv)")
    
    # 可选参数
    parser.add_argument('-o', '--out', default="merged_cami_input.tsv", help="[可选] 输出的合并文件名称 (默认: merged_cami_input.tsv)")
    parser.add_argument('--pct', action='store_true', help="[可选] 开启此标志后，会将 True_Abundance 列的数值乘以 100 转换为百分比格式。")

    args = parser.parse_args()
    
    merge_cami_data(args.gold, args.ref, args.out, args.pct)

if __name__ == "__main__":
    main()
