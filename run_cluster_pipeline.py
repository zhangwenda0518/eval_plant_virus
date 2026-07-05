#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import subprocess
import time
try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False
from collections import defaultdict

try:
    import polars as pl
except ImportError:
    sys.exit("❌ 致命错误: 找不到 polars 库。请先运行 `pip install polars` 安装。")

def parse_args():
    parser = argparse.ArgumentParser(description="🚀 终极聚类流水线: 自动运行聚类 / 现成结果提取 -> Polars极速分析 -> 提取代表与单簇")
    parser.add_argument("-i", "--fasta", required=True, help="输入的原始 FASTA 序列文件")
    parser.add_argument("-o", "--out_dir", default="Clustering_Pipeline_Out", help="输出根目录 (默认: Clustering_Pipeline_Out)")
    parser.add_argument("-m", "--method", choices=["mmseqs", "vclust", "both"], default="both", 
                        help="选择聚类工具类型。注意：当开启 --only_extract 时，必须明确指定为 'mmseqs' 或 'vclust'，不能为 'both'")
    
    # 聚类运行参数（非仅提取模式下生效）
    parser.add_argument("-t", "--threads", default="64", help="[聚类运行] 使用的线程数 (默认: 64)")
    parser.add_argument("--id", type=float, default=0.95, help="[聚类运行] 相似度/ANI 阈值 (默认: 0.95)")
    parser.add_argument("--cov", type=float, default=0.85, help="[聚类运行] 覆盖度阈值 (默认: 0.85)")
    parser.add_argument("--mmseqs-args", default="", help="[聚类运行] 透传给 MMseqs2 的额外参数 (空格分隔)")
    parser.add_argument("--vclust-algo", choices=['single','complete','uclust','cd-hit','set-cover','leiden'],
                        default="leiden", help="[聚类运行] VCLUST 聚类算法 (默认: leiden)")
    
    # 🎯 新增：仅提取与分析参数
    parser.add_argument("--only_extract", action="store_true", help="🔥 开启此选项则直接跳过聚类运行，仅对 --cluster_file 指定的现有结果进行分析与拆分")
    parser.add_argument("-c", "--cluster_file", default=None, help="现有的聚类结果文件路径 (当开启 --only_extract 时必填)")
    
    # 提取开关
    parser.add_argument("--skip_extract", action="store_true", help="跳过生成物理 FASTA 拆分文件 (默认不跳过，即默认执行释放提取)")
    
    return parser.parse_args()

def _tracked_run(tool, cmd, log_dir, shell=False):
    """执行子进程并记录时间/内存"""
    t0 = time.time()
    subprocess.run(cmd, shell=shell, check=True)
    t1 = time.time()
    wall_s = t1 - t0
    mem_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss if HAS_RESOURCE else 0
    os.makedirs(log_dir, exist_ok=True)
    log = os.path.join(log_dir, f"{tool}.time.mem.log")
    with open(log, 'a') as f:
        f.write(f"Time:{wall_s:.2f} seconds\nMemory:{mem_kb} KB\n")
    print(f"    ⏱ {wall_s:.1f}s, {mem_kb/1024:.0f} MB")
    return wall_s

def load_fasta_info(fasta_path):
    """极速加载 FASTA 序列并记录长度信息"""
    print(f"⏳ [阶段 1] 正在将 FASTA 文件读入内存字典...")
    fasta_dict = {}
    current_id, current_seq, current_header = None, [], ""
    
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_id:
                    seq_str = "".join(current_seq)
                    fasta_dict[current_id] = {"header": current_header, "seq": seq_str, "length": len(seq_str)}
                current_header = line
                current_id = line[1:].split()[0].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seq_str = "".join(current_seq)
            fasta_dict[current_id] = {"header": current_header, "seq": seq_str, "length": len(seq_str)}
            
    print(f"   ✅ 成功加载 {len(fasta_dict):,} 条序列。")
    return fasta_dict

def compute_and_export_stats_polars(clusters, fasta_dict, engine_out_dir):
    """利用 Polars 极速统计并输出 TSV 表格"""
    print(f"   -> 🚀 启动 Polars 聚合引擎计算单簇指标...")
    basic_stats = []

    for cname, info in clusters.items():
        ref_id = info["ref"]
        members = info["members"]

        ref_len = fasta_dict.get(ref_id, {}).get("length", 0)
        other_members = [m for m in members if m != ref_id]
        other_lens = [fasta_dict.get(m, {}).get("length", 0) for m in other_members if m in fasta_dict]

        basic_stats.append({
            "Cluster_ID": cname,
            "Ref_ID": ref_id,
            "Ref_Length": ref_len,
            "Total_Size": len(members),
            "Other_Members_Count": len(other_members),
            "Other_Max_Len": max(other_lens) if other_lens else 0,
            "Other_Avg_Len": round(sum(other_lens) / len(other_lens), 2) if other_lens else 0.0,
            "Other_Min_Len": min(other_lens) if other_lens else 0
        })

    schema_overrides = {"Other_Max_Len": pl.Int64, "Other_Min_Len": pl.Int64, "Other_Avg_Len": pl.Float64}
    df_final = pl.DataFrame(basic_stats, schema_overrides=schema_overrides)

    # 🎯 修复排序：提取纯数字作为临时列(sort_key)进行排序，排完再删掉临时列
    df_final = df_final.with_columns(
        pl.col("Cluster_ID").str.extract(r"(\d+)").cast(pl.Int64).alias("sort_key")
    ).sort("sort_key").drop("sort_key")

    report_file = os.path.join(engine_out_dir, "cluster_summary.tsv")
    df_final.write_csv(report_file, separator='\t')
    print(f"   ✅ 单簇统计报告已生成: {report_file}")

    return df_final

def generate_global_summary(df_final, total_input_seqs, engine_out_dir, engine_name):
    """提取数据框信息，输出全局概览并打印"""
    total_clusters = df_final.height
    cluster_sizes = df_final.get_column("Total_Size")
    
    max_size = cluster_sizes.max() if total_clusters > 0 else 0
    avg_size = cluster_sizes.mean() if total_clusters > 0 else 0
    singletons = df_final.filter(pl.col("Total_Size") == 1).height
    multi_clusters = total_clusters - singletons
    
    summary_text = (
        "===========================================================\n"
        f"📈 {engine_name.upper()} 全局聚类统计概览\n"
        "===========================================================\n"
        f"  📌 输入序列总数 (Total Input):  {total_input_seqs:,}\n"
        f"  📊 总聚类簇数 (Total Clusters): {total_clusters:,}\n"
        f"  📉 整体去冗余率 (Compression):  {(1 - total_clusters/total_input_seqs)*100:.2f}%\n"
        "  ---------------------------------------------------------\n"
        f"  👤 单例簇数 (Singletons):       {singletons:,} (占总簇数 {singletons/total_clusters*100:.2f}%)\n"
        f"  👨‍👩‍👧‍👦 多序列簇数 (Multi-member):  {multi_clusters:,}\n"
        f"  👑 最大簇包含序列数 (Max Size): {max_size:,}\n"
        f"  📏 平均簇大小 (Avg Size):       {avg_size:.2f}\n"
        "===========================================================\n"
    )
    
    print("\n" + summary_text)
    with open(os.path.join(engine_out_dir, "global_summary.txt"), 'w', encoding='utf-8') as f:
        f.write(summary_text)

def extract_fasta_files(clusters, fasta_dict, engine_out_dir):
    """物理提取：总代表序列 + 各个簇独立序列"""
    print(f"   -> 🚀 正在执行物理提取，生成最终 FASTA 文件...")
    extract_dir = os.path.join(engine_out_dir, "Split_Fastas")
    os.makedirs(extract_dir, exist_ok=True)
    
    global_ref_path = os.path.join(engine_out_dir, "all.cluster.ref.fasta")
    
    with open(global_ref_path, 'w', encoding='utf-8') as f_global:
        for cname, info in clusters.items():
            ref_id = info["ref"]
            members = info["members"]
            
            # 1. 写入全局代表序列合集
            if ref_id in fasta_dict:
                f_global.write(f"{fasta_dict[ref_id]['header']}\n{fasta_dict[ref_id]['seq']}\n")
                
            # 2. 写入独立的单簇代表序列
            with open(os.path.join(extract_dir, f"{cname}.ref.fasta"), 'w', encoding='utf-8') as fref:
                if ref_id in fasta_dict:
                    fref.write(f"{fasta_dict[ref_id]['header']}\n{fasta_dict[ref_id]['seq']}\n")
                    
            # 3. 写入独立的单簇所有序列
            with open(os.path.join(extract_dir, f"{cname}.all.fasta"), 'w', encoding='utf-8') as fall:
                for mid in members:
                    if mid in fasta_dict:
                        fall.write(f"{fasta_dict[mid]['header']}\n{fasta_dict[mid]['seq']}\n")

    print(f"   ✅ 提取完成！文件存入: {engine_out_dir}/")

def process_results(clusters, fasta_dict, engine_name, base_out_dir, skip_extract):
    """统一分析与提取调配器"""
    engine_out_dir = os.path.join(base_out_dir, f"{engine_name}_results")
    os.makedirs(engine_out_dir, exist_ok=True)
    print(f"\n⏳ [阶段 3 - {engine_name.upper()}] 开始剖析与提取结果...")
    
    # 计算统计并生成表
    df_final = compute_and_export_stats_polars(clusters, fasta_dict, engine_out_dir)
    # 生成全局概览
    generate_global_summary(df_final, len(fasta_dict), engine_out_dir, engine_name)
    # 物理提取
    if not skip_extract:
        extract_fasta_files(clusters, fasta_dict, engine_out_dir)
    else:
        print(f"   ℹ️ 接收到 --skip_extract，跳过文件物理拆分环节。")

def parse_existing_cluster_file(cluster_file, method):
    """🎯 用于 --only_extract 模式下直接解析现有的聚类结果文件"""
    print(f"⏳ [阶段 2 - 仅提取模式] 正在直接解析现有的 {method.upper()} 聚类文件...")
    raw_clusters = defaultdict(list)
    
    with open(cluster_file, 'r', encoding='utf-8') as f:
        if method == "vclust":
            f.readline()  # 跳过 Vclust 表头 'object \t cluster'
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    raw_clusters[parts[1].strip()].append(parts[0].strip())
            
            formatted_clusters = {}
            for cid, members in raw_clusters.items():
                formatted_clusters[f"cluster_{cid}"] = {"ref": members[0], "members": list(set(members))}
        else:  # mmseqs
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    raw_clusters[parts[0].strip()].append(parts[1].strip())
            
            formatted_clusters = {}
            for idx, (ref, members) in enumerate(raw_clusters.items()):
                formatted_clusters[f"cluster_{idx}"] = {"ref": ref, "members": list(set([ref] + members))}
                
    return formatted_clusters

def run_mmseqs(fasta, base_out_dir, threads, seq_id, cov, fasta_dict, skip_extract, mmseqs_args=""):
    print("\n" + "="*65)
    extra_str = f", Args: {mmseqs_args}" if mmseqs_args else ""
    print(f"🚀 [引擎] 启动 MMseqs2 easy-linclust 聚类 (ID: {seq_id}, Cov: {cov}{extra_str})")
    print("="*65)

    tmp_run_dir = os.path.join(base_out_dir, "tmp_mmseqs_run")
    os.makedirs(tmp_run_dir, exist_ok=True)
    out_prefix = os.path.join(tmp_run_dir, "mmseqs")

    cmd = [
        "mmseqs", "easy-linclust", fasta, out_prefix, os.path.join(tmp_run_dir, "tmp"),
        "--min-seq-id", str(seq_id), "--seq-id-mode", "1",
        "-c", str(cov), "-v", "0", "--threads", str(threads)
    ]
    # 默认参数（允许 --mmseqs-args 覆盖）
    extra_defaults = ["--cluster-mode", "2", "--cov-mode", "1"]
    if mmseqs_args:
        user_args = mmseqs_args.split()
        # 去除用户指定了的默认键
        keys = set(user_args[::2])
        filtered = []
        for i in range(0, len(extra_defaults), 2):
            if extra_defaults[i] not in keys:
                filtered.extend(extra_defaults[i:i+2])
        cmd.extend(filtered)
        cmd.extend(user_args)
    else:
        cmd.extend(extra_defaults)
    _tracked_run("mmseqs", cmd, tmp_run_dir)

    cluster_tsv = f"{out_prefix}_cluster.tsv"
    raw_clusters = defaultdict(list)
    if os.path.exists(cluster_tsv):
        with open(cluster_tsv, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    raw_clusters[parts[0].strip()].append(parts[1].strip())
                    
    formatted_clusters = {}
    for idx, (ref, members) in enumerate(raw_clusters.items()):
        formatted_clusters[f"cluster_{idx}"] = {"ref": ref, "members": list(set([ref] + members))}
        
    process_results(formatted_clusters, fasta_dict, "mmseqs", base_out_dir, skip_extract)

def run_vclust(fasta, base_out_dir, threads, ani, qcov, fasta_dict, skip_extract, algo="leiden"):
    print("\n" + "="*65)
    print(f"🚀 [引擎] 启动 Vclust 聚类 (ANI: {ani}, QCov: {qcov}, Algo: {algo})")
    print("="*65)
    
    tmp_run_dir = os.path.join(base_out_dir, "tmp_vclust_run")
    os.makedirs(tmp_run_dir, exist_ok=True)
    
    fltr_txt = os.path.join(tmp_run_dir, "fltr.txt")
    ani_tsv = os.path.join(tmp_run_dir, "ani.tsv")
    ani_ids = os.path.join(tmp_run_dir, "ani.ids.tsv")
    clusters_tsv = os.path.join(tmp_run_dir, "clusters.tsv")
    
    _tracked_run("vclust", f"vclust prefilter -i {fasta} -o {fltr_txt} --min-ident {ani} --threads {threads} -v 0", tmp_run_dir, shell=True)
    _tracked_run("vclust", f"vclust align -i {fasta} -o {ani_tsv} --filter {fltr_txt} --threads {threads} -v 0 ", tmp_run_dir, shell=True)
    _tracked_run("vclust", f"vclust cluster -i {ani_tsv} -o {clusters_tsv} --ids {ani_ids} --algorithm {algo} --metric ani --ani {ani} --qcov {qcov} -v 0 ", tmp_run_dir, shell=True)
    
    raw_clusters = defaultdict(list)
    if os.path.exists(clusters_tsv):
        with open(clusters_tsv, 'r') as f:
            f.readline()
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    raw_clusters[parts[1].strip()].append(parts[0].strip())
                    
    formatted_clusters = {}
    for cid, members in raw_clusters.items():
        formatted_clusters[f"cluster_{cid}"] = {"ref": members[0], "members": list(set(members))}
        
    process_results(formatted_clusters, fasta_dict, "vclust", base_out_dir, skip_extract)

def main():
    args = parse_args()
    start_time = time.time()
    
    # 基础安全性检查
    if args.only_extract:
        if not args.cluster_file:
            sys.exit("❌ 错误: 当开启 --only_extract 时，必须通过 -c / --cluster_file 指定现有的聚类文件！")
        if not os.path.exists(args.cluster_file):
            sys.exit(f"❌ 错误: 找不到指定的聚类文件: {args.cluster_file}")
        if args.method == "both":
            sys.exit("❌ 错误: 当开启 --only_extract 提取现有文件时，--method 必须明确为 'mmseqs' 或 'vclust'，不能为 'both'")
            
    os.makedirs(args.out_dir, exist_ok=True)
    
    print("======================================================================")
    print(f"🧬 聚类全栈流水线启动 | 模式: {'[仅提取与分析]' if args.only_extract else args.method.upper()}")
    print("==================================================================")
    
    # 阶段 1：统一加载 FASTA 序列信息
    fasta_dict = load_fasta_info(args.fasta)
    
    # 阶段 2：执行逻辑分流
    if args.only_extract:
        # 🎯 仅提取分支：直接解析传入的现有文件
        formatted_clusters = parse_existing_cluster_file(args.cluster_file, args.method)
        process_results(formatted_clusters, fasta_dict, args.method, args.out_dir, args.skip_extract)
    else:
        # 正常全自动运行分支
        if args.method in ["mmseqs", "both"]:
            run_mmseqs(args.fasta, args.out_dir, args.threads, args.id, args.cov, fasta_dict, args.skip_extract, args.mmseqs_args)
            
        if args.method in ["vclust", "both"]:
            run_vclust(args.fasta, args.out_dir, args.threads, args.id, args.cov, fasta_dict, args.skip_extract, args.vclust_algo)
        
    print("\n" + "="*66)
    print(f"🎉 任务圆满完成！总耗时: {time.time() - start_time:.2f} 秒")
    print(f"📁 结果已完好保存至: {os.path.abspath(args.out_dir)}")
    print("="*66)

if __name__ == "__main__":
    main()
