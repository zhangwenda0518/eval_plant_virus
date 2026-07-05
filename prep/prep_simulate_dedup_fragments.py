#!/usr/bin/env python3
"""
去冗余评估 — 模拟片段生成（非均匀分布）
参数：
  - 突变率：0%, 5%
  - 片段数量按长度比例非均匀分配：
    80%, 90% → 3-5 条
    60%, 70% → 10-20 条
    40%, 50% → 30-40 条
    20%, 30% → 50-60 条
  - >2500bp 的片段自动滑动窗口拆分为 2500bp 重叠子片段

每病毒生成约 440 条基础片段，拆分后约 600-800 条序列
500 病毒 × ~700 = ~350,000 条序列，FASTA 约 400 MB

用法:
  python prep_simulate_dedup_fragments.py \
      -i step1_dedup_viruses/ \
      -o step2_dedup_fragments.fasta --seed 42
"""

import argparse, os, random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── 参数配置 ──
MUTATION_RATES = [0.0, 0.05]
# 长度比例 → (最小片段数, 最大片段数)
FRAGMENT_QUOTA = [
    (0.9, 3, 5),
    (0.8, 3, 5),
    (0.7, 10, 20),
    (0.6, 10, 20),
    (0.5, 30, 40),
    (0.4, 30, 40),
    (0.3, 50, 60),
    (0.2, 50, 60),
]
WINDOW_SIZE = 2500
WINDOW_STEP = 1250


def introduce_mutations(seq_str, mutation_rate, rng):
    if mutation_rate == 0.0:
        return seq_str
    bases = ['A', 'T', 'C', 'G']
    mutated = list(seq_str)
    for i in range(len(mutated)):
        if rng.random() < mutation_rate:
            current = mutated[i].upper()
            choices = [b for b in bases if b != current]
            if choices:
                mutated[i] = rng.choice(choices)
    return "".join(mutated)


def apply_sliding_window(seq, base_id):
    """>2500bp 拆分，≤2500bp 原样返回"""
    seq_len = len(seq)
    if seq_len <= WINDOW_SIZE:
        return [SeqRecord(Seq(seq), id=base_id, description="")]
    records = []
    for i in range(0, seq_len - WINDOW_SIZE + 1, WINDOW_STEP):
        chunk = seq[i:i + WINDOW_SIZE]
        records.append(SeqRecord(Seq(chunk), id=f"{base_id}_sw{i}-{i+WINDOW_SIZE}", description=""))
    if (seq_len - WINDOW_SIZE) % WINDOW_STEP != 0:
        tail_start = seq_len - WINDOW_SIZE
        records.append(SeqRecord(Seq(seq[tail_start:]), id=f"{base_id}_sw{tail_start}-{seq_len}", description=""))
    return records


def simulate(record, rng, disable_split=False):
    """对一条参考基因组生成全部突变+片段化模拟序列"""
    genome_len = len(record.seq)
    species_id = record.id.split()[0]
    all_frags = []

    for mut_rate in MUTATION_RATES:
        mutated_full = introduce_mutations(str(record.seq), mut_rate, rng)
        mut_label = f"mut{int(mut_rate * 100)}pct"

        for frac, n_min, n_max in FRAGMENT_QUOTA:
            frag_len = int(genome_len * frac)
            if frag_len < 200:
                continue
            n_frags = rng.randint(n_min, n_max)
            frac_label = f"len{int(frac * 100)}pct"

            for i in range(n_frags):
                start = rng.randint(0, max(1, genome_len - frag_len))
                end = start + frag_len
                frag_seq = mutated_full[start:end]
                base_id = f"{species_id}_{mut_label}_{frac_label}_f{i+1}_pos{start}-{end}"
                if disable_split or len(frag_seq) <= WINDOW_SIZE:
                    all_frags.append(SeqRecord(Seq(frag_seq), id=base_id, description=""))
                else:
                    all_frags.extend(apply_sliding_window(frag_seq, base_id))

    return all_frags


def generate_host_fragments(host_genome, n_regions, host_label, rng):
    """从宿主基因组随机选取 n_regions 个 6k-20k 区间，再按 20%-90% 随机截断"""
    host_seqs = [(rec.id.split()[0], str(rec.seq)) for rec in SeqIO.parse(host_genome, "fasta")
                 if len(rec.seq) >= 20000]
    if not host_seqs:
        print("[host] WARNING: no sequences >= 20kb in host genome")
        return []

    host_frags = []
    for i in range(n_regions):
        _, seq = rng.choice(host_seqs)
        region_len = rng.randint(6000, min(20000, len(seq)))
        start = rng.randint(0, max(1, len(seq) - region_len))
        region = seq[start:start + region_len]
        frac = rng.uniform(0.2, 0.9)
        frag_len = int(region_len * frac)
        s = rng.randint(0, max(1, region_len - frag_len))
        frag = region[s:s + frag_len]
        base_id = f"{host_label}_host_f{i+1}_len{int(frac*100)}pct_pos{s}-{s+frag_len}"
        host_frags.extend(apply_sliding_window(frag, base_id))
    print(f"  [host] {n_regions} regions → {len(host_frags)} fragments")
    return host_frags


def main():
    parser = argparse.ArgumentParser(description="去冗余评估 — 模拟片段生成")
    parser.add_argument("-i", "--input", required=True,
                        help="选取的病毒 FASTA 目录或单个 FASTA 文件")
    parser.add_argument("-o", "--output", required=True,
                        help="输出的模拟 contig FASTA 文件")
    parser.add_argument("--host-genome", default=None,
                        help="宿主参考基因组 FASTA（生成负样本）")
    parser.add_argument("--n-host-regions", type=int, default=20,
                        help="宿主负样本区域数 (default: 20)")
    parser.add_argument("--host-label", default="HOST",
                        help="宿主片段 ID 前缀 (default: HOST)")
    parser.add_argument("--disable-split", action="store_true",
                        help="禁用 >2500bp 片段的滑动窗口拆分")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    if os.path.isdir(args.input):
        files = [os.path.join(args.input, f) for f in sorted(os.listdir(args.input))
                 if f.endswith(('.fasta', '.fa', '.fna'))]
    else:
        files = [args.input]

    all_frags = []
    for fpath in files:
        recs = list(SeqIO.parse(fpath, "fasta"))
        if not recs:
            continue
        rec = recs[0]  # 每文件一个基因组
        print(f"  {rec.id.split()[0]} ({len(rec.seq)} bp)")
        frags = simulate(rec, rng, args.disable_split)
        all_frags.extend(frags)
        print(f"    → {len(frags)} fragments")

    # 负样本：宿主基因组随机区域
    if args.host_genome and os.path.exists(args.host_genome):
        print(f"\n[host] Generating negative fragments from {args.host_genome}...")
        host_frags = generate_host_fragments(args.host_genome, args.n_host_regions,
                                             args.host_label, rng)
        all_frags.extend(host_frags)

    SeqIO.write(all_frags, args.output, "fasta")
    n_species = len(set(f.split('_')[0] for f in os.listdir(args.input) if f.endswith('.fasta')))
    if args.host_genome:
        n_species += 1
    print(f"\n✅ 共生成 {len(all_frags)} 条序列 ({n_species} species + host) → {args.output}")
    print(f"   文件大小: {os.path.getsize(args.output) / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
