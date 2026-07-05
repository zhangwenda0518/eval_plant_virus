#!/bin/bash
# ============================================================
# 评估六：序列去重聚类 — 统一运行脚本
# 金标准: 每条片段 ID 中的物种名 (如 NC_002030_mut0pct_len50pct_f1 → NC_002030)
#
# 用法: bash run_dedup_clustering.sh --input step2_dedup_fragments.fasta --outdir step3_dedup_cluster/
# ============================================================
set -e

INPUT="step2_dedup_fragments.fasta"
OUTDIR="step3_dedup_cluster"
THREADS=30
MIN_ID=0.90

show_help() {
    echo "Usage: bash run_dedup_clustering.sh [options]"
    echo "  --input FILE     模拟片段 FASTA"
    echo "  --outdir DIR     输出目录"
    echo "  --threads N      线程数 (default: 30)"
    echo "  --min-id FLOAT   CD-HIT/VCLUST 相似度阈值 (default: 0.90)"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --input)   INPUT="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --min-id)  MIN_ID="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) echo "Unknown: $1"; show_help ;;
    esac
done

mkdir -p "$OUTDIR"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

N_SEQS=$(grep -c '^>' "$INPUT")
log "Input: $INPUT ($N_SEQS sequences)"

# ── 1. MMseqs2 easy-cluster ──
log "=== MMseqs2 ==="
MMSEQS_OUT="$OUTDIR/mmseqs2"
if [ ! -f "$MMSEQS_OUT/.DONE" ]; then
    mkdir -p "$MMSEQS_OUT"
    mmseqs easy-cluster "$INPUT" "$MMSEQS_OUT/mmseqs_result" "$MMSEQS_OUT/tmp" \
        --min-seq-id "$MIN_ID" -c 0.8 --cov-mode 0 --threads "$THREADS" \
        --cluster-mode 2 2>&1 | tail -5
    touch "$MMSEQS_OUT/.DONE"
    log "  MMseqs2 done"
else
    log "  MMseqs2 skip (.DONE)"
fi

# ── 2. CD-HIT ──
log "=== CD-HIT ==="
CDHIT_OUT="$OUTDIR/cdhit"
if [ ! -f "$CDHIT_OUT/.DONE" ]; then
    mkdir -p "$CDHIT_OUT"
    cd-hit-est -i "$INPUT" -o "$CDHIT_OUT/cdhit_result" \
        -c "$MIN_ID" -n 10 -d 0 -M 32000 -T "$THREADS" 2>&1 | tail -5
    touch "$CDHIT_OUT/.DONE"
    log "  CD-HIT done"
else
    log "  CD-HIT skip (.DONE)"
fi

# ── 3. VCLUST (vsearch --cluster_fast) ──
log "=== VCLUST ==="
VCLUST_OUT="$OUTDIR/vclust"
if [ ! -f "$VCLUST_OUT/.DONE" ]; then
    mkdir -p "$VCLUST_OUT"
    vsearch --cluster_fast "$INPUT" \
        --id "$MIN_ID" --centroids "$VCLUST_OUT/vclust_centroids.fasta" \
        --uc "$VCLUST_OUT/vclust_result.uc" --threads "$THREADS" 2>&1 | tail -5
    touch "$VCLUST_OUT/.DONE"
    log "  VCLUST done"
else
    log "  VCLUST skip (.DONE)"
fi

# ── 4. dRep (依赖 checkm, 对核苷酸聚类) ──
# dRep 需要每个 genome 独立文件，这里用 --cluster_fast 替代
# 如果 dRep 不可用，用 mmseqs linclust 替代
log "=== dRep / MMseqs linclust ==="
DREP_OUT="$OUTDIR/linclust"
if [ ! -f "$DREP_OUT/.DONE" ]; then
    mkdir -p "$DREP_OUT"
    mmseqs easy-linclust "$INPUT" "$DREP_OUT/linclust_result" "$DREP_OUT/tmp" \
        --min-seq-id "$MIN_ID" -c 0.8 --cov-mode 0 --threads "$THREADS" 2>&1 | tail -5
    touch "$DREP_OUT/.DONE"
    log "  linclust done"
else
    log "  linclust skip (.DONE)"
fi

log "=== All clustering done ==="
log "Output: $OUTDIR"
