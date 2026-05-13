#!/usr/bin/env python3
import os
import glob
import re
import subprocess
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 解决 SSH 环境下无 GUI 导致的报错
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==========================================
# 1. 全局配置与样式定义
# ==========================================
MY_PALETTE = {
    "Megahit": "#4C72B0",
    "Penguin": "#DD8452",
    "RNAViralSPAdes": "#55A868",
    "RefineC_Merge": "#C44E52"
}
MUT_ORDER = ["0.0pct", "5.0pct", "15.0pct", "30.0pct"]
TOOL_ORDER = ["Megahit", "Penguin", "RNAViralSPAdes", "RefineC_Merge"]

def map_tool_name(name_str):
    """统一不同文件命名中的工具名称"""
    name_lower = str(name_str).lower()
    if "megahit" in name_lower: return "Megahit"
    elif "penguin" in name_lower: return "Penguin"
    elif "rnaviralspades" in name_lower: return "RNAViralSPAdes"
    elif "refinec" in name_lower or "merged" in name_lower: return "RefineC_Merge"
    return name_str

def get_existing_file(base_path):
    """智能匹配 .fasta 或 .fasta.gz"""
    for p in [base_path, base_path + ".gz"]:
        if os.path.exists(p): return p
    return None

# ==========================================
# 2. 评估执行模块 (MetaQUAST)
# ==========================================
def run_single_metaquast(target_dir, args):
    """单个样本的 MetaQUAST 任务"""
    path_parts = os.path.normpath(target_dir).split(os.sep)
    dataset_folder = path_parts[0] if len(path_parts) > 1 else "."
    base_name = path_parts[-1]
    
    out_dir = os.path.join(dataset_folder, args.out_dir, base_name)
    
    # 待评估的文件基础路径
    expected = [
        os.path.join(target_dir, f"{base_name}_megahit.contig.fasta"),
        os.path.join(target_dir, f"{base_name}_penguin.contig.fasta"),
        os.path.join(target_dir, f"{base_name}_rnaviralspades.contig.fasta"),
        os.path.join(target_dir, f"{base_name}_all_tools_refineC_merge.merged.fasta")
    ]

    actual_files = []
    for f in expected:
        path = get_existing_file(f)
        if path: actual_files.append(path)

    if len(actual_files) == 0:
        return f"❌ 失败: {base_name} (未找到任何组装文件)"

    cmd = [
        "metaquast", "-o", out_dir, "-r", args.ref_dir,
        "--min-contig", str(args.min_contig),
        "-t", str(args.threads),
        "-l", "Megahit,Penguin,RNAViralSPAdes,RefineC_Merge"
    ] + actual_files

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return f"✅ 成功: [{dataset_folder}] {base_name}"
    except Exception:
        return f"❌ 运行出错: {base_name}"

# ==========================================
# 3. 数据解析模块
# ==========================================
def parse_all_data(args):
    print("\n📥 正在提取 MetaQUAST 质量指标...")
    # 扫描 Dataset_Mut_*/metaquast_results/*/runs_per_reference/*/transposed_report.tsv
    report_pattern = os.path.join(args.input_pattern.split('/')[0], args.out_dir, "*", "runs_per_reference", "*", "transposed_report.tsv")
    report_files = glob.glob(report_pattern)
    
    q_data = []
    for path in report_files:
        try:
            parts = path.split(os.sep)
            mut_rate = parts[0].replace("Dataset_Mut_", "")
            sample_name = parts[2]
            virus_name = parts[4]
            depth = int(sample_name.split("_")[3]) # 假设 Master_Sample_Depth_100...
            
            df = pd.read_csv(path, sep='\t')
            df['Tool'] = df['Assembly'].apply(map_tool_name)
            df['Mutation_Rate'] = mut_rate
            df['Reads_Subsample'] = depth
            df['Virus'] = virus_name
            q_data.append(df)
        except: continue
    
    df_quality = pd.concat(q_data, ignore_index=True) if q_data else pd.DataFrame()

    print("📥 正在提取资源消耗日志 (.time.mem.log)...")
    log_files = glob.glob(os.path.join(args.input_pattern, "*.time.mem.log"))
    res_data = []
    for path in log_files:
        try:
            parts = path.split(os.sep)
            with open(path, 'r') as f:
                content = f.read()
                t_match = re.search(r'Time:\s*([\d\.]+)', content, re.I)
                m_match = re.search(r'Memory:\s*(\d+)', content, re.I)
                if t_match and m_match:
                    res_data.append({
                        "Mutation_Rate": parts[0].replace("Dataset_Mut_", ""),
                        "Reads_Subsample": int(parts[2].split("_")[3]),
                        "Tool": map_tool_name(parts[3]),
                        "Time_Seconds": float(t_match.group(1)),
                        "Memory_MB": float(m_match.group(1)) / 1024.0
                    })
        except: continue
    
    df_resource = pd.DataFrame(res_data)
    return df_quality, df_resource

# ==========================================
# 4. 绘图模块 (严格复现原 6 张图)
# ==========================================
def plot_benchmark(df_q, df_r):
    if df_q.empty: return
    sns.set_theme(style="ticks", font_scale=1.2)
    
    # 质量指标循环绘图
    metrics = [
        ("Genome fraction (%)", "Plot_01_Genome_Fraction.png", "Genome Fraction (%)"),
        ("Largest alignment", "Plot_02_Continuity.png", "Length (bp)"),
        ("# misassemblies", "Plot_03_Misassemblies.png", "Number of Errors"),
        ("# mismatches per 100 kbp", "Plot_04_Mismatches.png", "Mismatches / 100 kbp"),
        ("# indels per 100 kbp", "Plot_05_Indels.png", "Indels / 100 kbp")
    ]

    for col, fname, ylabel in metrics:
        if col not in df_q.columns: continue
        g = sns.catplot(
            data=df_q, x="Reads_Subsample", y=col, hue="Tool",
            row="Mutation_Rate", col="Virus", kind="box",
            row_order=[m for m in MUT_ORDER if m in df_q['Mutation_Rate'].unique()],
            hue_order=TOOL_ORDER, palette=MY_PALETTE, height=4, aspect=1.3
        )
        g.fig.suptitle(f"Benchmark: {col}", fontsize=20, y=1.02, fontweight='bold')
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   📊 已保存: {fname}")

    if not df_r.empty:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sns.lineplot(data=df_r, x="Reads_Subsample", y="Time_Seconds", hue="Tool", marker="o", palette=MY_PALETTE, ax=axes[0])
        axes[0].set_title("A. Runtime (Seconds)", fontweight='bold')
        sns.lineplot(data=df_r, x="Reads_Subsample", y="Memory_MB", hue="Tool", marker="s", palette=MY_PALETTE, ax=axes[1])
        axes[1].set_title("B. Peak Memory (MB)", fontweight='bold')
        plt.tight_layout()
        plt.savefig("Plot_06_Resource_Consumption.png", dpi=300)
        plt.close()
        print(f"   📊 已保存: Plot_06_Resource_Consumption.png")

# ==========================================
# 5. 主程序控制
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="🧬 病毒组装 Benchmark 一体化流水线")
    parser.add_argument("-i", "--input_pattern", default="Dataset_Mut_*/1.virus-assembly/Master_*", help="输入目录模式")
    parser.add_argument("-r", "--ref_dir", default="../virus-db", help="参考库目录")
    parser.add_argument("-o", "--out_dir", default="metaquast_results", help="MetaQUAST 输出子目录名")
    parser.add_argument("-j", "--jobs", type=int, default=4, help="并行样本数")
    parser.add_argument("-t", "--threads", type=int, default=4, help="每个任务线程数")
    parser.add_argument("-m", "--min_contig", type=int, default=100)
    parser.add_argument("--skip-run", action="store_true", help="跳过执行，仅分析绘图")
    
    args = parser.parse_args()

    # 第一阶段：执行 MetaQUAST
    if not args.skip_run:
        dirs = glob.glob(args.input_pattern)
        if not dirs:
            print("❌ 错误: 未找到符合条件的目录！")
            return
        print(f"🚀 开始批量评估 {len(dirs)} 个样本...")
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(run_single_metaquast, d, args): d for d in dirs}
            for f in as_completed(futures): print(f.result())

    # 第二阶段：数据汇总与分析
    df_q, df_r = parse_all_data(args)
    
    if not df_q.empty:
        df_q.to_csv("benchmark_quality_summary.csv", index=False)
        plot_benchmark(df_q, df_r)
        print("\n🎉 全部任务执行完毕！")
    else:
        print("\n⚠️ 未提取到有效数据，请检查路径。")

if __name__ == "__main__":
    main()
