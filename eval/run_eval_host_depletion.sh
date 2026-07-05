#!/bin/bash
# ============================================================
# 评估一：宿主过滤消融实验 — 完整调度脚本
# 支持断点续跑、时间/内存记录
#
# 用法:
#   bash run_eval_host_depletion.sh
#   bash run_eval_host_depletion.sh --threads 20 --jobs 20
# ============================================================
set -e
BIN_DIR="$(cd "$(dirname "$0")/.." && pwd)/deps"

# ---- 默认值 ----
LOGDIR="step5_host_free_logs"
HOST_K2="$HOME/database/host_db/kraken2"
HOST_HISAT2="$HOME/database/host_db/hisat2/host"
INPUT_DIR="step5_LoD_all"
OUTDIR="step5_host_free"
SEQ_TYPE="rna-short"
TOOL="hisat2"
JOBS=10
THREADS=10
FILTER="true"
SKIP_D0=false
SKIP_D1=false
SKIP_D2=false
SKIP_D3=false
SKIP_D4=false

# ---- 解析命名参数 ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --logdir)        LOGDIR="$2"; shift 2 ;;
        --host-k2)       HOST_K2="$2"; shift 2 ;;
        --host-hisat2)   HOST_HISAT2="$2"; shift 2 ;;
        --input-dir)     INPUT_DIR="$2"; shift 2 ;;
        --output-dir)    OUTDIR="$2"; shift 2 ;;
        --seq-type)      SEQ_TYPE="$2"; shift 2 ;;
        --tool)          TOOL="$2"; shift 2 ;;
        --jobs)          JOBS="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --filter)        FILTER="$2"; shift 2 ;;
        --skip-d0)       SKIP_D0=true; shift ;;
        --skip-d1)       SKIP_D1=true; shift ;;
        --skip-d2)       SKIP_D2=true; shift ;;
        --skip-d3)       SKIP_D3=true; shift ;;
        --skip-d4)       SKIP_D4=true; shift ;;
        -h|--help)
            echo "Usage: bash run_eval_host_depletion.sh [options]"
            echo ""
            echo "  宿主过滤消融实验: D0(基线) D1(仅Kraken2) D2(仅HISAT2) D3(K2+HISAT2) D4(完整三步)"
            echo ""
            echo "Options:"
            echo "  --host-k2 DIR      Kraken2宿主数据库"
            echo "  --host-hisat2 DIR  HISAT2宿主索引前缀"
            echo "  --input-dir DIR    LoD数据目录"
            echo "  --output-dir DIR   输出根目录 (default: step5_host_free)"
            echo "  --seq-type STR     测序类型 (default: rna-short)"
            echo "  --jobs N           并发样本数 (default: 10)"
            echo "  --threads N        单样本线程 (default: 10)"
            echo "  --skip-d0 ~ d4     跳过指定实验组"
            exit 0
            ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

mkdir -p "$LOGDIR" "$OUTDIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/pipeline.log"
}

run_with_log() {
    local name="$1"; shift
    local outdir="$1"; shift
    local logfile="$LOGDIR/${name}.log"

    if [ -f "$outdir/.DONE" ]; then
        log "[$name] Already completed, skipping"
        return 0
    fi

    log "[$name] Starting... (log: $logfile)"
    local t0; t0=$(date +%s)
    /usr/bin/time -v -o "$LOGDIR/${name}.time" "$@" > "$logfile" 2>&1
    local ret=$?
    local t1; t1=$(date +%s)
    local elapsed=$((t1 - t0))
    local mem_peak
    mem_peak=$(grep "Maximum resident set size" "$LOGDIR/${name}.time" 2>/dev/null | awk '{print $NF}')

    if [ $ret -eq 0 ]; then
        touch "$outdir/.DONE"
        log "[$name] Done | Time: ${elapsed}s | Peak RSS: ${mem_peak:-N/A} KB"
    else
        log "[$name] Failed (exit code: $ret) | Time: ${elapsed}s"
        return $ret
    fi
}

# ---- D0: 基线（无过滤）----
if [ "$SKIP_D0" != true ]; then
    D0_DIR="${OUTDIR}/D0_baseline"
    if [ ! -f "$D0_DIR/.DONE" ]; then
        log "[D0] Creating baseline (symlink raw LoD data)..."
        mkdir -p "$D0_DIR"
        for f in "$INPUT_DIR"/*_R1.fastq.gz; do
            base=$(basename "$f" _R1.fastq.gz)
            [ ! -e "$D0_DIR/${base}_R1.fastq.gz" ] && ln -sf "$(realpath "$f")" "$D0_DIR/"
            r2="${f/_R1/_R2}"
            [ -e "$r2" ] && [ ! -e "$D0_DIR/${base}_R2.fastq.gz" ] && ln -sf "$(realpath "$r2")" "$D0_DIR/"
        done
        touch "$D0_DIR/.DONE"
    fi
    log "[D0] Baseline: $(ls "$D0_DIR"/*_R1.fastq.gz 2>/dev/null | wc -l) samples"
fi

# ---- D1: 仅 Kraken2 ----
if [ "$SKIP_D1" != true ]; then
    run_with_log "D1_kraken2_only" "${OUTDIR}/D1_kraken2_only" \
        python "$BIN_DIR/host_depletion.py" \
            --seq-type "$SEQ_TYPE" --tool "$TOOL" \
            -k "$HOST_K2" -x "$HOST_HISAT2" \
            -I "$INPUT_DIR" -O "${OUTDIR}/D1_kraken2_only/" \
            --steps kraken2 --jobs "$JOBS" --threads "$THREADS" \
            --tmp "/tmp/host_free_D1_$$" -f "$FILTER" \
            --logs_dir "${LOGDIR}/D1_kraken2_only/" || true
fi

# ---- D2: 仅 HISAT2 ----
if [ "$SKIP_D2" != true ]; then
    run_with_log "D2_hisat2_only" "${OUTDIR}/D2_hisat2_only" \
        python "$BIN_DIR/host_depletion.py" \
            --seq-type "$SEQ_TYPE" --tool "$TOOL" \
            -k "$HOST_K2" -x "$HOST_HISAT2" \
            -I "$INPUT_DIR" -O "${OUTDIR}/D2_hisat2_only/" \
            --steps align --jobs "$JOBS" --threads "$THREADS" \
            --tmp "/tmp/host_free_D2_$$" -f "$FILTER" \
            --logs_dir "${LOGDIR}/D2_hisat2_only/" || true
fi

# ---- D3: Kraken2 + HISAT2 ----
if [ "$SKIP_D3" != true ]; then
    run_with_log "D3_k2_hisat2" "${OUTDIR}/D3_k2_hisat2" \
        python "$BIN_DIR/host_depletion.py" \
            --seq-type "$SEQ_TYPE" --tool "$TOOL" \
            -k "$HOST_K2" -x "$HOST_HISAT2" \
            -I "$INPUT_DIR" -O "${OUTDIR}/D3_k2_hisat2/" \
            --steps kraken2,align --jobs "$JOBS" --threads "$THREADS" \
            --tmp "/tmp/host_free_D3_$$" -f "$FILTER" \
            --logs_dir "${LOGDIR}/D3_k2_hisat2/" || true
fi

# ---- D4: 完整三步 (Kraken2 + HISAT2 + rRNA) ----
if [ "$SKIP_D4" != true ]; then
    run_with_log "D4_full" "${OUTDIR}/D4_full" \
        python "$BIN_DIR/host_depletion.py" \
            --seq-type "$SEQ_TYPE" --tool "$TOOL" \
            -k "$HOST_K2" -x "$HOST_HISAT2" \
            -I "$INPUT_DIR" -O "${OUTDIR}/D4_full/" \
            --steps kraken2,align,rrna --rrna \
            --jobs "$JOBS" --threads "$THREADS" \
            --tmp "/tmp/host_free_D4_$$" -f "$FILTER" \
            --logs_dir "${LOGDIR}/D4_full/" || true
fi

# ---- 汇总 ----
log "====================="
log "Host Depletion Ablation Summary"
log "====================="
for cfg in D0_baseline D1_kraken2_only D2_hisat2_only D3_k2_hisat2 D4_full; do
    time_file="$LOGDIR/${cfg}.time"
    done_mark="${OUTDIR}/${cfg}/.DONE"
    if [ -f "$done_mark" ]; then
        if [ -f "$time_file" ]; then
            elapsed=$(grep "Elapsed" "$time_file" 2>/dev/null | awk '{print $NF}' || echo "N/A")
            rss=$(grep "Maximum resident set size" "$time_file" 2>/dev/null | awk '{print $NF}' || echo "N/A")
            echo "✅ $cfg | Time=$elapsed | Peak RSS=${rss}KB"
        else
            echo "✅ $cfg"
        fi
    else
        echo "❌ $cfg | Not completed"
    fi
done | tee -a "$LOGDIR/pipeline.log"

log "Pipeline finished."
