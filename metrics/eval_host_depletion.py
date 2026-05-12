#!/usr/bin/env python3
"""
eval_host_depletion.py (终极自动检测版)
— 自动探测 D0~D4 目录、无强制时间限制、多进程并行。
"""

import os, sys, re, glob, argparse, subprocess
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from tqdm import tqdm

def parse_filename(basename):
    """
    智能提取 Accession 和 LoD_Factor，兼容 _PE_R1 和 _PE_clean_R1
    """
    m = re.search(r'LoD_Mixed_(.*?)_(\d+\.\d+)_', basename)
    if m:
        return m.group(1), float(m.group(2))
    return None, None


def count_reads(fastq_file, virus_acc):
    """
    精确统计：完全取消超时限制，依靠底层系统命令本身的报错机制。
    """
    virus_count, total_count = 0, 0
    basename = os.path.basename(fastq_file)

    # 1. 统计病毒 Reads
    try:
        vc = subprocess.run(
            ["rg", "-z", "-c", f"^@{virus_acc}-", fastq_file],
            capture_output=True, text=True
        )
        if vc.stdout.strip().isdigit():
            virus_count = int(vc.stdout.strip())
        elif vc.stderr:
            print(f"\n[Error] rg 读取 {basename} 时报错: {vc.stderr.strip()}")
    except Exception as e:
        print(f"\n[Error] 扫描病毒 Reads 失败 {basename}: {e}")

    # 2. 统计总 Reads
    try:
        tc_cmd = f"zcat '{fastq_file}' 2>/dev/null | wc -l"
        tc = subprocess.run(tc_cmd, shell=True, capture_output=True, text=True)

        if not tc.stdout.strip().isdigit() or int(tc.stdout.strip()) == 0:
            tc_cmd = f"gzcat '{fastq_file}' 2>/dev/null | wc -l"
            tc = subprocess.run(tc_cmd, shell=True, capture_output=True, text=True)

        if tc.stdout.strip().isdigit():
            total_count = int(tc.stdout.strip()) // 4

        if tc.stderr:
            print(f"\n[Error] zcat 读取 {basename} 时报错: {tc.stderr.strip()}")

    except Exception as e:
        print(f"\n[Error] 统计总 Reads 失败 {basename}: {e}")

    if total_count == 0:
        print(f"\n[Alert] 发现 0 Reads 文件 (文件可能为空或损坏): {basename}")

    host_count = max(0, total_count - virus_count)
    return host_count, virus_count


def process_file(args_tuple):
    fastq_file, group_name = args_tuple
    basename = os.path.basename(fastq_file)
    virus_acc, lod_factor = parse_filename(basename)

    if not virus_acc:
        return None

    host, virus = count_reads(fastq_file, virus_acc)
    total = host + virus

    return {
        'Group': group_name,
        'File': basename,
        'Accession': virus_acc,
        'LoD_Factor': lod_factor,
        'Host_Reads': host,
        'Virus_Reads': virus,
        'Total_Reads': total,
        'Virus_Pct': round(virus / total * 100, 4) if total > 0 else 0,
    }


def scan_dir(dir_path, group_name, jobs):
    files = sorted(glob.glob(os.path.join(dir_path, "*_R1.fastq*")))
    tasks = [(f, group_name) for f in files]

    if not tasks:
        print(f"\n[Warning] 在 {dir_path} 中没有找到任何 FASTQ 文件！")
        return []

    print(f"[{group_name}] 正在使用 {jobs} 个并发进程扫描 {len(tasks)} 个文件...")
    records = []

    with ProcessPoolExecutor(max_workers=jobs) as executor:
        results = list(tqdm(executor.map(process_file, tasks), total=len(tasks), desc=group_name))

    for res in results:
        if res is not None:
            records.append(res)

    total_host = sum(r['Host_Reads'] for r in records)
    total_virus = sum(r['Virus_Reads'] for r in records)
    print(f"  -> 完成! 提取到 {len(records)} 个样本记录. 总宿主: {total_host:,} | 总病毒: {total_virus:,}\n")
    return records


def main():
    parser = argparse.ArgumentParser(description="宿主过滤消融实验评估 (自动输入探测版)")
    # 将繁琐的 --d0~d4 替换为一个核心输入目录
    parser.add_argument('-i', '--input_dir', required=True, help='包含 D0~D4 文件夹的父目录 (如: step5_host_free)')
    parser.add_argument('--outdir', default='step5_host_free_analysis', help='输出目录')
    parser.add_argument('-j', '--jobs', type=int, default=8, help='并行运行的进程数')
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"[Error] 输入目录不存在: {args.input_dir}")
        sys.exit(1)

    # ================= 自动探测逻辑 =================
    dir_map = {'D0': None, 'D1': None, 'D2': None, 'D3': None, 'D4': None}

    for folder in os.listdir(args.input_dir):
        folder_path = os.path.join(args.input_dir, folder)
        if not os.path.isdir(folder_path):
            continue

        # 截取文件夹前两个字符并转大写，用于匹配 D0~D4
        prefix = folder[:2].upper()
        if prefix in dir_map:
            dir_map[prefix] = folder_path

    # 校验是否找齐了所有必须的文件夹
    missing = [k for k, v in dir_map.items() if v is None]
    if missing:
        print(f"[Error] 探测失败！在 {args.input_dir} 中找不到以下前缀的文件夹: {', '.join(missing)}")
        print("请确保文件夹以 D0_, D1_, D2_, D3_, D4_ 开头。")
        sys.exit(1)

    print("\n✅ 成功检测到目录结构:")
    for k in sorted(dir_map.keys()):
        print(f"  {k} -> {os.path.basename(dir_map[k])}")
    print("-" * 50 + "\n")
    # ================================================

    os.makedirs(args.outdir, exist_ok=True)

    all_records = []
    # 这里的 group_name (如 'D0_Baseline') 是为了保证和后续的绘图脚本强绑定，不能改。
    all_records.extend(scan_dir(dir_map['D0'], 'D0_Baseline', args.jobs))
    all_records.extend(scan_dir(dir_map['D1'], 'D1_Kraken2_only', args.jobs))
    all_records.extend(scan_dir(dir_map['D2'], 'D2_HISAT2_only', args.jobs))
    all_records.extend(scan_dir(dir_map['D3'], 'D3_K2_HISAT2', args.jobs))
    all_records.extend(scan_dir(dir_map['D4'], 'D4_Full', args.jobs))

    detail_df = pd.DataFrame(all_records)

    d0_df = detail_df[detail_df['Group'] == 'D0_Baseline']
    d0_by_key = {(r['Accession'], r['LoD_Factor']): r['Virus_Reads'] for _, r in d0_df.iterrows()}

    def calc_retention(row):
        if row['Group'] == 'D0_Baseline': return 100.0
        base = d0_by_key.get((row['Accession'], row['LoD_Factor']), None)
        return round(row['Virus_Reads'] / base * 100, 2) if base and base > 0 else None

    detail_df['Virus_Retention'] = detail_df.apply(calc_retention, axis=1)

    detail_path = os.path.join(args.outdir, 'host_depletion_detail.tsv')
    detail_df.to_csv(detail_path, sep='\t', index=False)

    print(f"🎉 数据处理完毕！详细结果已保存至: {detail_path}")


if __name__ == "__main__":
    main()