#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import polars as pl
import argparse
import os
import sys
import time

# ==========================================
# 常量字典
# ==========================================
DOMAIN_KEYWORDS = {
    "Bacteria": ["bacteria"],
    "Archaea": ["archaea"],
    "Fungi":["fungi", "fungus"],
    "Viruses":["virus", "viruses", "viria", "virae", "viricota", "viricetes", "virales", "viridae", "virinae"]
}

RANK_DICT = {
    "P": "phylum", "C": "class", "O": "order", 
    "F": "family", "G": "genus", "S": "species"
}

# ==========================================
# 1. 树节点类
# ==========================================
class TreeNode:
    def __init__(self, taxid, p_taxid, rank, level_num, name, track_tools):
        self.taxid = taxid
        self.p_taxid = p_taxid
        self.rank = rank
        self.level_num = level_num
        self.name = name
        self.parent = None
        self.children = []
        
        self.lvl_reads = {t: 0 for t in track_tools} # L3 (Direct)
        self.all_reads = {t: 0 for t in track_tools} # L2 (Clade)
        
        self.lineage = ""
        self.species = ""
        self.genus = ""
        self.domain = ""

    def add_child(self, node):
        self.children.append(node)

# ==========================================
# 2. 建树与特征推导
# ==========================================
def build_taxonomy_tree(tax_file, track_tools):
    taxid2node = {}
    root_nodes = []
    
    print(f"[1/6] Loading taxonomy tree from {tax_file}...")
    with open(tax_file, 'r') as f:
        for line in f:
            if not line.strip(): continue
            parts = [x.strip() for x in line.split('|')]
            if len(parts) >= 5:
                try:
                    tid, pid, rank, lvl_num, name = int(parts[0]), int(parts[1]), parts[2], int(parts[3]), parts[4]
                    taxid2node[tid] = TreeNode(tid, pid, rank, lvl_num, name, track_tools)
                except ValueError:
                    continue
                    
    for tid, node in taxid2node.items():
        if tid == node.p_taxid or node.p_taxid == 0 or tid == 1:
            root_nodes.append(node)
        elif node.p_taxid in taxid2node:
            p_node = taxid2node[node.p_taxid]
            node.parent = p_node
            p_node.add_child(node)
            
    queue = root_nodes.copy()
    while queue:
        curr = queue.pop(0)
        
        if curr.parent:
            curr.lineage = curr.parent.lineage + "|" + curr.name
            curr.genus = curr.name if curr.rank in ['G', 'genus'] else curr.parent.genus
            curr.species = curr.name if curr.rank in ['S', 'species'] else curr.parent.species
        else:
            curr.lineage = curr.name
            curr.genus = curr.name if curr.rank in ['G', 'genus'] else ""
            curr.species = curr.name if curr.rank in ['S', 'species'] else ""
            
        assigned_dom = "unclassified"
        name_lower = curr.name.lower()
        for dom_key, kw_list in DOMAIN_KEYWORDS.items():
            if any(kw in name_lower for kw in kw_list):
                assigned_dom = dom_key
                break
                
        if assigned_dom != "unclassified":
            curr.domain = assigned_dom
        elif curr.parent and curr.parent.domain:
            curr.domain = curr.parent.domain
        else:
            curr.domain = "unclassified"
            
        for child in curr.children:
            queue.append(child)
            
    for node in taxid2node.values():
        if not node.species: node.species = node.name

    return taxid2node, root_nodes

# ==========================================
# 3. 数据加载与合并算法 (极致内存优化版)
# ==========================================
def load_tool_output(file_path, tool_id):
    # 构建空惰性表，占位符 (防止文件缺失时报错)
    # 【优化】使用 UInt32 存储 TaxID，将数字矩阵的内存消耗减少 50%
    empty_lf = pl.LazyFrame(schema={"Seq_ID": pl.String, f"{tool_id}_tax_ID": pl.UInt32})

    # 检查文件是否不存在或为空 (0 bytes)
    if not file_path or not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        print(f"  [Warning] {tool_id} output '{file_path}' is empty or missing. Skipping.")
        return empty_lf

    if tool_id == 'centrifuger': seq_col, tax_col = "column_1", "column_3"
    elif tool_id == 'metabuli': seq_col, tax_col = "column_2", "column_3"
    else: seq_col, tax_col = "column_2", "column_3"
        
    try:
        # 【优化】使用 scan_csv 惰性读取，不立即载入内存
        lf = pl.scan_csv(file_path, separator='\t', has_header=False, ignore_errors=True, truncate_ragged_lines=True)
        
        lf = lf.select([
            pl.col(seq_col).cast(pl.String).str.replace(r"/[12]$", "").alias("Seq_ID"),
            pl.col(tax_col).cast(pl.UInt32, strict=False).alias(f"{tool_id}_tax_ID")
        ]).filter(pl.col(f"{tool_id}_tax_ID").is_not_null())
        
        return lf
    except Exception as e:
        # 捕获异常，比如列名不匹配、格式损坏，保证脚本继续运行
        print(f"  [Warning] Failed to scan {tool_id} output '{file_path}': {e}. Skipping.")
        return empty_lf

def get_lca_for_taxids(valid_tids, taxid2node):
    if not valid_tids: return 0
    if len(valid_tids) == 1: return valid_tids[0]
    current_lca = valid_tids[0]
    for tid in valid_tids[1:]:
        ancestors = set()
        curr = current_lca
        while curr != 0 and curr != 1 and curr in taxid2node:
            ancestors.add(curr)
            curr = taxid2node[curr].p_taxid
        ancestors.add(1) 
        curr = tid
        while curr not in ancestors:
            if curr not in taxid2node or curr == 0:
                curr = 1; break
            curr = taxid2node[curr].p_taxid
        current_lca = curr
    return current_lca

def get_lowest_or_lca(tids, taxid2node):
    valid_tids = list(set([t for t in tids if t != 0 and t in taxid2node]))
    if not valid_tids: return 0
    if len(valid_tids) == 1: return valid_tids[0]
    
    valid_tids.sort(key=lambda t: taxid2node[t].level_num, reverse=True)
    deepest_tid = valid_tids[0]
    
    is_linear = True
    for other_tid in valid_tids[1:]:
        curr = deepest_tid
        found = False
        while curr != 0 and curr != 1 and curr in taxid2node:
            if curr == other_tid: found = True; break
            curr = taxid2node[curr].p_taxid
        if not found:
            is_linear = False; break
            
    return deepest_tid if is_linear else get_lca_for_taxids(valid_tids, taxid2node)

# ==========================================
# 主流程
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Metagenomics Consensus Builder & Report Analyzer")
    parser.add_argument('--kraken2', help="Raw Kraken2 .out file")
    parser.add_argument('--kaiju', help="Raw Kaiju .out file")
    parser.add_argument('--metabuli', help="Raw Metabuli tsv file")
    parser.add_argument('--centrifuger', help="Raw Centrifuger classed.tsv file")
    parser.add_argument('--kunpeng', help="Raw Kun-peng output file")
    parser.add_argument('--krakenuniq', help="Raw KrakenUniq .txt/.out file")
    parser.add_argument('--kraken2x', help="Raw Kraken2-x .out file")
    
    parser.add_argument('-tax', required=True, help="Taxonomy file from make_ktaxonomy.py")
    parser.add_argument('-outdir', default=".", help="Output directory")
    parser.add_argument('-prefix', default="sample", help="Output prefix")
    
    parser.add_argument('--mode', choices=['lowest', 'lca', 'prefer', 'all'], default='lowest', 
                        help="Merging strategy: 'lowest', 'lca', 'prefer', or 'all' (skip consensus).")
    parser.add_argument('--prefer', type=str, help="Required if --mode is 'prefer' (e.g., kraken2)")
    
    args = parser.parse_args()

    tools_provided = []
    tool_files = {}
    SUPPORTED_TOOLS = ['kraken2', 'kaiju', 'metabuli', 'centrifuger', 'kunpeng', 'krakenuniq', 'kraken2x']
    for tool in SUPPORTED_TOOLS:
        if getattr(args, tool):
            tools_provided.append(tool)
            tool_files[tool] = getattr(args, tool)

    if not tools_provided:
        print("ERROR: Provide at least one tool output.")
        return
    if args.mode == 'prefer' and (not args.prefer or args.prefer not in tools_provided):
        print(f"ERROR: You must specify a valid --prefer tool when using '--mode prefer'. Available: {tools_provided}")
        return

    if not os.path.exists(args.outdir): os.makedirs(args.outdir)
    start_time = time.time()

    track_tools = tools_provided.copy() if args.mode == 'all' else tools_provided + ["Consensus"]
    taxid2node, root_nodes = build_taxonomy_tree(args.tax, track_tools)

    print("[2/6] Parsing output files & building Read Matrix (Lazy & Streaming)...")
    vrl_lazy = None
    for tool in tools_provided:
        lf = load_tool_output(tool_files[tool], tool)
        if vrl_lazy is None:
            vrl_lazy = lf
        else:
            # 增量构造连接树，coalesce=True 避免 Seq_ID 出现多列
            vrl_lazy = vrl_lazy.join(lf, on="Seq_ID", how="full", coalesce=True)

    # 填充空值为 0 (unclassified)，并全部统一为 UInt32 以极致节省内存
    fill_exprs = [pl.col(f"{t}_tax_ID").fill_null(0).cast(pl.UInt32) for t in tools_provided]
    vrl_lazy = vrl_lazy.with_columns(fill_exprs)

    # 【核心：触发计算】
    print("      -> Collecting lazy frames into memory. Streaming to prevent OOM...")
    try:
        # 开启 streaming=True，分块合并防止内存溢出
        vrl = vrl_lazy.collect(streaming=True)
    except Exception as e:
        print(f"      -> [Warning] Streaming failed ({e}), falling back to standard collect...")
        vrl = vrl_lazy.collect()

    # 增加对极端情况的保护：所有工具的数据全是空的
    if vrl.height == 0:
        print(f"\n[Warning] All tool outputs for {args.prefix} are empty. Generating blank reports and exiting gracefully.")
        with open(os.path.join(args.outdir, f'{args.prefix}_tree_report.tsv'), 'w') as f:
            f.write("No mapped reads found.\n")
        sys.exit(0) # 平稳退出，不影响批量处理的下一个样本

    # ===============================================
    # 根据模式决断 (Arbitration)
    # ===============================================
    if args.mode != 'all':
        print(f"[3/6] Applying Consensus Strategy: '{args.mode.upper()}' ...")
        tool_cols = [f"{t}_tax_ID" for t in tools_provided]
        
        if args.mode == 'prefer':
            vrl = vrl.with_columns(pl.col(f"{args.prefer}_tax_ID").alias("Consensus_tax_ID"))
        else:
            # 提取唯一组合计算 LCA / Lowest，极大加快速度
            unique_patterns = vrl.select(tool_cols).unique().to_dicts()
            for row in unique_patterns:
                tids = [row[c] for c in tool_cols]
                if args.mode == 'lowest':
                    row['Consensus_tax_ID'] = get_lowest_or_lca(tids, taxid2node)
                elif args.mode == 'lca':
                    valid_tids = [t for t in tids if t != 0]
                    row['Consensus_tax_ID'] = get_lca_for_taxids(valid_tids, taxid2node)
                    
            # Schema 显式指定 Consensus_tax_ID 为 UInt32
            consensus_schema = {c: pl.UInt32 for c in tool_cols}
            consensus_schema['Consensus_tax_ID'] = pl.UInt32
            
            consensus_df = pl.DataFrame(unique_patterns, schema=consensus_schema)
            vrl = vrl.join(consensus_df, on=tool_cols, how="left")
    else:
        print("[3/6] Mode 'ALL' selected: Skipping read-level consensus arbitration...")

    print("[4/6] Aggregating multi-tool counts up the tree...")
    def ensure_node(tid):
        if tid not in taxid2node:
            name = "unclassified" if tid == 0 else f"TaxID_{tid}"
            node = TreeNode(tid, 0, 'U', 0, name, track_tools)
            node.domain, node.lineage = "unclassified", name
            node.species, node.genus = name, ""
            taxid2node[tid] = node
            if tid != 0: root_nodes.append(node)

    for tool in track_tools:
        tool_cnt = vrl.filter(pl.col(f"{tool}_tax_ID") != 0).group_by(f"{tool}_tax_ID").agg(pl.len().alias("count")).to_dicts()
        for row in tool_cnt:
            tid, count = row[f"{tool}_tax_ID"], row["count"]
            ensure_node(tid)
            node = taxid2node[tid]
            node.lvl_reads[tool] += count
            node.all_reads[tool] += count
            p = node.parent
            while p:
                p.all_reads[tool] += count
                p = p.parent

    totals = {t: len(vrl) for t in track_tools}
    unclassified_counts = {t: vrl.filter(pl.col(f"{t}_tax_ID") == 0).height for t in track_tools}

    print("[5/6] Generating standard .report files and Side-by-Side Tree View...")
    for tool in track_tools:
        rep_path = os.path.join(args.outdir, f"{args.prefix}_{tool}.report")
        with open(rep_path, 'w') as f:
            if unclassified_counts[tool] > 0:
                pct = (unclassified_counts[tool] / totals[tool]) * 100 if totals[tool] else 0.0
                f.write(f"{pct:6.2f}\t{unclassified_counts[tool]}\t{unclassified_counts[tool]}\tU\t0\tunclassified\n")
            def dfs_rep(node):
                if node.all_reads[tool] == 0: return
                pct = (node.all_reads[tool] / totals[tool]) * 100 if totals[tool] else 0.0
                f.write(f"{pct:6.2f}\t{node.all_reads[tool]}\t{node.lvl_reads[tool]}\t{node.rank}\t{node.taxid}\t{'  '*node.level_num}{node.name}\n")
                for child in sorted(node.children, key=lambda c: c.all_reads[tool], reverse=True): dfs_rep(child)
            for r in sorted(root_nodes, key=lambda c: c.all_reads[tool], reverse=True): dfs_rep(r)

    tree_rows = []
    def get_sort_key(node):
        return node.all_reads["Consensus"] if args.mode != 'all' else sum(node.all_reads[t] for t in tools_provided)
        
    def build_tree_out(node):
        if all(node.all_reads[t] == 0 for t in track_tools): return
        row = {}
        for t in track_tools:
            row[f"Pct_{t}"] = (node.all_reads[t] / totals[t] * 100) if totals[t] else 0.0
        for t in track_tools:
            row[f"Clade_{t}"] = node.all_reads[t]
        for t in track_tools:
            row[f"Direct_{t}"] = node.lvl_reads[t]
            
        row["Rank"] = node.rank
        row["TaxID"] = node.taxid
        row["Indented_Name"] = "  " * node.level_num + node.name
        
        tree_rows.append(row)
        for child in sorted(node.children, key=get_sort_key, reverse=True):
            build_tree_out(child)

    if 0 in taxid2node: build_tree_out(taxid2node[0])
    for r in sorted(root_nodes, key=get_sort_key, reverse=True): build_tree_out(r)
    
    if tree_rows:
        pl.DataFrame(tree_rows).with_columns([pl.col(f"Pct_{t}").round(2) for t in track_tools]) \
                               .write_csv(os.path.join(args.outdir, f'{args.prefix}_tree_report.tsv'), separator='\t')

    print("[6/6] Generating flat Rank & Lineage Summary...")
    master_records = []
    for node in taxid2node.values():
        if all(node.all_reads[t] == 0 for t in track_tools): continue
        rec = {"domain": node.domain, "rank": node.rank, "tax_id": node.taxid, "name": node.name, "lineage": node.lineage}
        for t in track_tools:
            rec[f"reads_clade_{t}"] = node.all_reads[t]
            rec[f"reads_direct_{t}"] = node.lvl_reads[t]
        master_records.append(rec)
    
    if master_records:
        df_master = pl.DataFrame(master_records)

        def write_multidim_csv(domain_val, rank_val, out_file):
            df_sub = df_master.filter(pl.col("domain") == domain_val) if domain_val else df_master
            df_sub = df_sub.filter(pl.col("rank") == rank_val) if rank_val else df_sub
            if df_sub.height == 0: return
                
            rel_exprs = []
            for t in track_tools:
                t_tot = df_sub[f"reads_clade_{t}"].sum()
                expr = pl.when(t_tot > 0).then(pl.col(f"reads_clade_{t}") / t_tot * 100).otherwise(0.0).round(4).alias(f"relative_abundance_{t}")
                rel_exprs.append(expr)
                
            df_sub = df_sub.with_columns(rel_exprs)
            
            cols = ["domain", "rank", "tax_id", "name"]
            if args.mode != 'all':
                cols.extend(["relative_abundance_Consensus", "reads_clade_Consensus", "reads_direct_Consensus"])
            for t in tools_provided:
                cols.extend([f"relative_abundance_{t}", f"reads_clade_{t}", f"reads_direct_{t}"])
            cols.append("lineage")
            
            if args.mode != 'all':
                df_sub = df_sub.select(cols).sort(by="reads_clade_Consensus", descending=True)
            else:
                sort_cols = [f"reads_clade_{t}" for t in tools_provided]
                df_sub = df_sub.select(cols).sort(by=sort_cols, descending=[True]*len(sort_cols))
                
            df_sub.write_csv(out_file)

        for dom in DOMAIN_KEYWORDS.keys():
            for r_code, r_name in RANK_DICT.items():
                write_multidim_csv(dom, r_code, os.path.join(args.outdir, f"{args.prefix}_{dom.lower()}_{r_name}.csv"))
        write_multidim_csv(None, "S", os.path.join(args.outdir, f"{args.prefix}_all_species.csv"))

    vrl.write_csv(os.path.join(args.outdir, f'{args.prefix}_VRL_Mapping.txt'), separator='\t')

    print(f"\n[Success] Processed {args.prefix} ({args.mode.upper()}) in {time.time() - start_time:.2f} seconds!")

if __name__ == "__main__":
    main()
