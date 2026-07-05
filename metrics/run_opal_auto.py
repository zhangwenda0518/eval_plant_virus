#!/usr/bin/env python3
import os
import glob
import subprocess
import sys
import argparse

def main():
    # === 设置命令行参数解析 ===
    parser = argparse.ArgumentParser(
        description="🚀 CAMI 数据自动化合并与 OPAL 评估工具\n自动提取、合并指定目录下的 CAMI profile 并一键运行 OPAL 评估。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('-i', '--input-dir', default="ALL.cami_profiles", 
                        help="包含各个软件预测结果(.profile)的输入目录\n(默认: ALL.cami_profiles)")
    parser.add_argument('-g', '--gold-dir', default="eval_cami_gold_results", 
                        help="包含金标准(.profile)的输入目录\n(默认: eval_cami_gold_results)")
    parser.add_argument('-o', '--out-dir', default="OPAL_Results_ALL_Datasets", 
                        help="OPAL 评估报告生成的输出目录\n(默认: OPAL_Results_ALL_Datasets)")
    parser.add_argument('--merged-gold', default="ALL_GoldStandard.profile", 
                        help="合并后的金标准中间文件名\n(默认: ALL_GoldStandard.profile)")
    
    args = parser.parse_args()

    # === 路径配置 (应用参数) ===
    gold_dir = args.gold_dir
    cami_dir = args.input_dir
    output_gold = args.merged_gold
    opal_out_dir = args.out_dir

    print(f"🚀 开始进行 CAMI 数据自动化合并与 OPAL 评估...")
    print(f"📁 金标准目录: {gold_dir}")
    print(f"📁 预测集目录: {cami_dir}")
    print(f"📁 报告输出到: {opal_out_dir}\n")

    # ==========================================
    # 步骤 1：合并所有的金标准 (Gold Standards)
    # ==========================================
    gold_files = glob.glob(f"{gold_dir}/*_GoldStandard.profile")
    if not gold_files:
        print(f"⚠️ 警告: 在 {gold_dir} 目录下未找到金标准文件！请检查路径。")
    else:
        print(f"📦 正在合并 {len(gold_files)} 个金标准文件 -> {output_gold} ...")
        with open(output_gold, "w", encoding="utf-8") as outfile:
            for fname in gold_files:
                with open(fname, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    outfile.write(content)
                    # 确保文件末尾有换行符，防止合并错行
                    if not content.endswith('\n'):
                        outfile.write('\n')

    # ==========================================
    # 步骤 2：动态解析软件列表并合并 Profiles
    # ==========================================
    profile_files = glob.glob(f"{cami_dir}/*.profile")
    if not profile_files:
        print(f"❌ 错误: 在 {cami_dir} 目录下找不到任何 .profile 文件！")
        sys.exit(1)

    tools = set()
    for fname in profile_files:
        basename = os.path.basename(fname)
        # 假设文件名格式为 <sample_name>_<tool>.profile
        name_without_ext = basename.replace('.profile', '')
        tool_name = name_without_ext.split('_')[-1]
        tools.add(tool_name)

    # 排序以保证每次运行顺序一致
    tools = sorted(list(tools))
    print(f"🔍 自动识别到 {len(tools)} 个评估工具: {', '.join(tools)}")

    combined_profiles = []
    for tool in tools:
        tool_files = glob.glob(f"{cami_dir}/*_{tool}.profile")
        out_profile = f"ALL_{tool}.profile"
        combined_profiles.append(out_profile)
        
        print(f"   => 合并 {tool} 的 {len(tool_files)} 个样本 -> {out_profile} ...")
        with open(out_profile, "w", encoding="utf-8") as outfile:
            for fname in tool_files:
                with open(fname, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    outfile.write(content)
                    if not content.endswith('\n'):
                        outfile.write('\n')

    # ==========================================
    # 步骤 3：运行 OPAL
    # ==========================================
    # 将工具首字母大写作为 OPAL 的 Label (例如 kraken2 -> Kraken2)
    labels = [tool.capitalize() for tool in tools]
    labels_str = ",".join(labels)

    print("\n📊 正在启动 OPAL 进行多维度评估评分...")
    
    # 构建 OPAL 命令
    opal_cmd = [
        "opal.py", 
        "-g", output_gold, 
        "-o", opal_out_dir, 
        "--labels", labels_str
    ] + combined_profiles

    # 打印即将运行的命令（方便调试）
    print(f"▶️  执行命令: {' '.join(opal_cmd)}\n")

    try:
        subprocess.run(opal_cmd, check=True)
        print(f"\n🎉 完美！评估报告已生成在 {opal_out_dir} 目录下！")
    except subprocess.CalledProcessError as e:
        print(f"\n❌ OPAL 运行失败，返回码: {e.returncode}")
    except FileNotFoundError:
        print("\n❌ 找不到 opal.py，请确保 OPAL 已正确安装并加入了环境变量！")

if __name__ == "__main__":
    main()
