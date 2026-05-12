#!/usr/bin/env python3
"""
提取含病毒保守结构域远缘同源的宿主蛋白编码序列（负样本B来源）

两阶段工作流：
  阶段1（自动）：调用 para_hmmscan.pl 对宿主蛋白进行 Pfam 注释
  阶段2（自动）：从注释结果中筛选含病毒保守结构域的序列，截取核酸片段

病毒保守结构域（Diamond BLASTx假阳性主要来源）:
  - DEAD/DEAH 解旋酶 (PF00270, PF00271)
  - 蛋白激酶 (PF00069)
  - 锌指蛋白 (PF00096, PF13912)
  - NTPase (PF00437)
  - RdRp远缘同源 (PF00680, PF00978)

用法:
  # 完整流程：蛋白FASTA → Pfam注释 → 筛选陷阱序列
  python prep_extract_conserved_traps.py \
      --host-proteins host_proteins.faa \
      --host-genome all.genome.uniq.fasta \
      --host-gff genomic.gff \
      --pfam-db ~/database/pfam-v35/Pfam-A.hmm \
      --outdir conserved_traps/ \
      --n-sequences 300

  # 如果已有 Pfam 注释结果，跳过阶段1
  python prep_extract_conserved_traps.py \
      --host-proteins host_proteins.faa \
      --host-genome all.genome.uniq.fasta \
      --host-gff genomic.gff \
      --pfam-annot host.pfam.tsv \
      --outdir conserved_traps/ \
      --n-sequences 300

输出:
  {outdir}/
    ├── pfam_annotation.tsv          # Pfam 注释结果
    ├── trap_0000.fasta ...          # 陷阱序列文件
    └── trap_summary.tsv             # 按 Pfam 类型统计
"""

import argparse, os, sys, random, csv, subprocess, shutil
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ============================================================
# 目标保守结构域 Pfam ID
# ============================================================
TARGET_PFAM = {
    "PF00270": "DEAD_helicase",
    "PF00271": "DEAH_helicase",
    "PF00069": "Protein_kinase",
    "PF00096": "Zinc_finger_C2H2",
    "PF13912": "Zinc_finger_C2H2_type",
    "PF00437": "NTPase_T2SS",
    "PF00680": "RdRp_1",
    "PF00978": "RdRp_2",
}


# ============================================================
# 阶段1: 运行 para_hmmscan.pl 进行 Pfam 注释
# ============================================================

def run_pfam_annotation(protein_fasta, pfam_db, out_dir, cpu=20, hmmscan_cpu=4, chunk=100):
    """
    调用 para_hmmscan.pl 对宿主蛋白序列进行 Pfam 注释。
    使用 --outformat 输出简化7列格式。
    """
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # 检查 pfam_db 是否已建索引
    for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
        if not os.path.exists(pfam_db + ext):
            print(f"[pfam] Building HMMER index for {pfam_db}...")
            subprocess.run(["hmmpress", pfam_db], check=True)
            break

    # 输出文件
    domtbl_out = os.path.join(out_dir, "pfam_annotation.tsv")

    if os.path.exists(domtbl_out) and os.path.getsize(domtbl_out) > 100:
        print(f"[pfam] Annotation already exists: {domtbl_out}")
        return domtbl_out

    # 构建 para_hmmscan.pl 命令
    cmd = [
        "perl", os.path.expanduser("~/bin/para_hmmscan.pl"),
        "--hmm_db", pfam_db,
        "--cpu", str(cpu),
        "--hmmscan_cpu", str(hmmscan_cpu),
        "--chunk", str(chunk),
        "--outformat",
        "--evalue1", "1e-5",
        "--evalue2", "1e-3",
        "--coverage", "0.25",
        "--tmp_prefix", os.path.join(out_dir, "hmmscan"),
        protein_fasta,
    ]

    print(f"[pfam] Running para_hmmscan.pl on {protein_fasta}...")
    print(f"[pfam] Command: {' '.join(cmd)}")

    with open(domtbl_out, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)

    print(f"[pfam] Pfam annotation saved to {domtbl_out}")
    return domtbl_out


# ============================================================
# 阶段2: 解析 Pfam 结果并提取陷阱序列
# ============================================================

def parse_pfam_output(pfam_tsv):
    """
    解析 para_hmmscan.pl --outformat 的简化7列输出:
      GeneID \t HMM_accession \t HMM_Name \t E-value \t Score \t Coverage \t Description

    也兼容标准 hmmscan --domtblout 格式。
    """
    pfam = defaultdict(list)
    with open(pfam_tsv, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            # 判断格式：简化7列 vs domtblout (23列)
            # Pfam ID 格式: PF开头后跟5位数字，可能带版本号.xx
            is_pfam_col = parts[1].startswith("PF") and parts[1][2:7].isdigit() if len(parts[1]) >= 7 else False
            if len(parts) >= 7 and is_pfam_col:
                # 简化7列格式: GeneID PFxxxxx HMM_Name E-value Score Coverage Description
                seq_id = parts[0].strip()
                pfam_id = parts[1].strip()
            else:
                # 尝试标准 domtblout 格式：第1列target name, 第4列Pfam accession
                seq_id = parts[0].strip()
                pfam_id = parts[3].strip() if len(parts) > 3 else ""
                if not pfam_id.startswith("PF"):
                    # 扫描全部列找PF ID
                    for cell in parts:
                        cell = cell.strip()
                        if cell.startswith("PF") and len(cell) == 7:
                            pfam_id = cell
                            break

            base_id = pfam_id.split(".")[0] if "." in pfam_id else pfam_id
            if base_id in TARGET_PFAM:
                pfam[seq_id].append(base_id)

    total = len(pfam)
    print(f"[parse] DEBUG: pfam dict has {total} keys, first 5: {list(pfam.keys())[:5]}")
    print(f"[parse] Found {total} sequences with target Pfam domains")
    if total > 0:
        for pid, desc in TARGET_PFAM.items():
            count = sum(1 for v in pfam.values() if pid in v)
            if count > 0:
                print(f"  {pid} ({desc}): {count}")
    return pfam


def get_cds_sequences(host_genome, host_gff, host_proteins, out_dir):
    """
    获取CDS核酸序列。优先级：
      1. 从 GFF3 + 基因组用 gffread 提取
      2. 直接加载宿主蛋白 FASTA（作为fallback）
    """
    cds_seqs = {}

    # 方法1: gffread
    if host_genome and host_gff and os.path.exists(host_gff):
        if shutil.which("gffread"):
            print("[cds] Extracting CDS from GFF with gffread...")
            cds_fa = os.path.join(out_dir, "extracted_cds.fasta")
            try:
                subprocess.run(
                    ["gffread", "-x", cds_fa, "-g", host_genome, host_gff],
                    check=True, capture_output=True)
                if os.path.exists(cds_fa):
                    for rec in SeqIO.parse(cds_fa, "fasta"):
                        cds_seqs[rec.id] = str(rec.seq)
                    print(f"[cds] Got {len(cds_seqs)} CDS from gffread")
            except Exception:
                print("[cds] gffread failed, falling back to protein FASTA")
        else:
            print("[cds] gffread not found, falling back to protein FASTA")

    # 方法2: 蛋白FASTA作为fallback
    # 蛋白ID通常对应转录本ID，去掉 .p1 等后缀后可能匹配CDS ID
    if host_proteins and os.path.exists(host_proteins):
        print("[cds] Loading protein sequences as fallback...")
        added = 0
        for rec in SeqIO.parse(host_proteins, "fasta"):
            if rec.id not in cds_seqs:
                cds_seqs[rec.id] = str(rec.seq)
                added += 1
        if added > 0:
            print(f"[cds] Added {added} protein sequences as CDS fallback")

    return cds_seqs


# ============================================================
# 主逻辑
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="提取保守结构域陷阱序列")
    parser.add_argument("--host-proteins", required=True,
                        help="宿主蛋白序列 FASTA")
    parser.add_argument("--host-genome", help="宿主基因组 FASTA")
    parser.add_argument("--host-gff", help="宿主基因组 GFF3 注释")
    parser.add_argument("--pfam-db", default=os.path.expanduser("~/database/pfam-v35/Pfam-A.hmm"),
                        help="Pfam-A.hmm 数据库路径")
    parser.add_argument("--pfam-annot", help="已有 Pfam 注释 TSV（跳过阶段1）")
    parser.add_argument("--cpu", type=int, default=12,
                        help="para_hmmscan.pl 并行数")
    parser.add_argument("--hmmscan-cpu", type=int, default=4,
                        help="单个 hmmscan 线程数")
    parser.add_argument("--chunk", type=int, default=100,
                        help="每块蛋白序列数")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--n-sequences", type=int, default=300)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # ---- 阶段1: Pfam 注释 ----
    if args.pfam_annot and os.path.exists(args.pfam_annot):
        print("[1/3] Using existing Pfam annotation...")
        pfam_tsv = args.pfam_annot
    else:
        print("[1/3] Running Pfam annotation via para_hmmscan.pl...")
        pfam_tsv = run_pfam_annotation(
            args.host_proteins, args.pfam_db, args.outdir,
            args.cpu, args.hmmscan_cpu, args.chunk)

    # ---- 阶段2: 解析 + 筛选 ----
    print("[2/3] Parsing Pfam results and filtering target domains...")
    target_seqs = parse_pfam_output(pfam_tsv)

    if not target_seqs:
        print("[WARNING] No sequences with target Pfam domains found. Creating placeholder.")
        for i in range(args.n_sequences):
            out_file = os.path.join(args.outdir, f"trap_{i:04d}.fasta")
            with open(out_file, "w") as f:
                f.write(f">trap_{i:04d}|source=synthetic|type=placeholder\n")
                f.write("".join(rng.choices("ACGT", k=rng.randint(800, 3000))) + "\n")
        print(f"[DONE] Created {args.n_sequences} placeholder sequences")
        return

    # ---- 阶段3: 提取核酸序列并写出 ----
    print("[3/3] Extracting nucleic acid sequences and writing output...")
    cds_seqs = get_cds_sequences(args.host_genome, args.host_gff,
                                  args.host_proteins, args.outdir)

    # 匹配 Pfam 命中的蛋白ID到 CDS 核酸序列
    qualified = []
    for seq_id, pfam_list in target_seqs.items():
        seq = cds_seqs.get(seq_id, "")
        if len(seq) >= 500:
            qualified.append((seq_id, seq, pfam_list))
        else:
            # 模糊匹配：蛋白ID可能带 .p1 后缀，或转换为CDS ID
            for k, v in cds_seqs.items():
                base_id = seq_id.split(".p")[0] if ".p" in seq_id else seq_id
                if base_id in k and len(v) >= 500:
                    qualified.append((k, v, pfam_list))
                    break

    # 去重（同一序列可能命中多个 Pfam）
    seen_ids = set()
    unique_qualified = []
    for q in qualified:
        if q[0] not in seen_ids:
            seen_ids.add(q[0])
            unique_qualified.append(q)

    print(f"[filter] {len(unique_qualified)} unique qualified sequences (>=500 bp)")

    if len(unique_qualified) == 0:
        print("[WARNING] No qualified sequences, creating placeholder")
        for i in range(args.n_sequences):
            out_file = os.path.join(args.outdir, f"trap_{i:04d}.fasta")
            with open(out_file, "w") as f:
                f.write(f">trap_{i:04d}|source=synthetic|type=placeholder\n")
                f.write("".join(rng.choices("ACGT", k=rng.randint(800, 3000))) + "\n")
        return

    # 按 Pfam 类型分层选取
    by_pfam = defaultdict(list)
    for seq_id, seq, pfam_list in unique_qualified:
        for pid in pfam_list:
            by_pfam[pid].append((seq_id, seq, pfam_list))

    selected = []
    if by_pfam:
        n_per_pfam = max(1, args.n_sequences // len(by_pfam))
        for pid, members in by_pfam.items():
            n_pick = min(n_per_pfam, len(members))
            picked = rng.sample(members, n_pick)
            selected.extend(picked)

    # 达到目标数量
    if len(selected) > args.n_sequences:
        selected = rng.sample(selected, args.n_sequences)
    elif len(selected) < args.n_sequences:
        # 不足时从剩余中补齐
        sel_set = {(s[0], s[1]) for s in selected}
        remaining = [(q[0], q[1], q[2]) for q in unique_qualified 
                     if (q[0], q[1]) not in sel_set]
        if remaining:
            extra = rng.sample(remaining,
                               min(args.n_sequences - len(selected), len(remaining)))
            for q in extra:
                selected.append(tuple(q))

    # 写出
    summary = defaultdict(int)
    for i, (seq_id, seq, pfam_list) in enumerate(selected):
        # 截取 500-3000 bp 片段（模拟组装的碎片化contig）
        frag_len = min(rng.randint(500, 3000), len(seq))
        start = rng.randint(0, max(1, len(seq) - frag_len))
        frag = seq[start:start + frag_len]

        out_file = os.path.join(args.outdir, f"trap_{i:04d}.fasta")
        with open(out_file, "w") as f:
            f.write(f">trap_{i:04d}|source={seq_id}|pfam={','.join(pfam_list)}\n{frag}\n")

        for pid in pfam_list:
            summary[TARGET_PFAM.get(pid, pid)] += 1

    # 写出摘要
    with open(os.path.join(args.outdir, "trap_summary.tsv"), "w") as f:
        f.write("pfam\tdescription\tcount\n")
        for desc, count in sorted(summary.items(), key=lambda x: -x[1]):
            f.write(f"{desc}\t{count}\n")

    print(f"[DONE] Created {len(selected)} trap sequences in {args.outdir}")
    print(f"  Pfam distribution:")
    for desc, count in sorted(summary.items(), key=lambda x: -x[1]):
        print(f"    {desc}: {count}")


if __name__ == "__main__":
    main()
