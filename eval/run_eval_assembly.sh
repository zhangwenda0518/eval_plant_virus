#!/bin/bash
# ============================================================
# 评估二：病毒组装方法比较 — 自动扫描 + 分批调度
# 自动识别模拟数据目录下的组和突变率
#
# 用法:
#   bash run_eval_assembly.sh
#   bash run_eval_assembly.sh --sim-data ./data --output-dir ./asm --batch 3 --threads 20 --jobs 5
# ============================================================
set -e

# ---- 默认值 ----
SIMDIR="step2_benchmark_data"
OUTDIR="step6_assemblies"
LOGDIR="step6_logs"
REF_FASTA="step6_ref_viruses.fasta"
VIRUS_DIR="step1_eval_viruses"
BATCH=2
THREADS=15
JOBS=4

show_help() {
    echo "Usage: bash run_eval_assembly.sh [options]"
    echo ""
    echo "  自动扫描 --sim-data 下的 group_*/Dataset_Mut_*pct/Jackknife_Subsamples/"
    echo "  分批运行 assembly_pipeline.py，避免内存溢出。.DONE 标记支持断点续跑。"
    echo ""
    echo "Options:"
    echo "  --sim-data DIR      模拟数据根目录 (default: step2_benchmark_data)"
    echo "  --output-dir DIR    组装输出目录 (default: step6_assemblies)"
    echo "  --logdir DIR        日志目录 (default: step6_logs)"
    echo "  --virus-dir DIR     病毒参考基因组目录 (default: step1_eval_viruses)"
    echo "  --ref-fasta FILE    合并后的参考FASTA (default: step6_ref_viruses.fasta)"
    echo "  --batch N           每批并发组数 (default: 2)"
    echo "  --threads N         单样本线程数 (default: 15)"
    echo "  --jobs N            assembly_pipeline 内部并发 (default: 4)"
    echo "  -h, --help          显示帮助"
    exit 0
}

# ---- 解析参数 ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --sim-data)   SIMDIR="$2"; shift 2 ;;
        --output-dir) OUTDIR="$2"; shift 2 ;;
        --logdir)     LOGDIR="$2"; shift 2 ;;
        --virus-dir)  VIRUS_DIR="$2"; shift 2 ;;
        --ref-fasta)  REF_FASTA="$2"; shift 2 ;;
        --batch)      BATCH="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --jobs)       JOBS="$2"; shift 2 ;;
        -h|--help)    show_help ;;
        *) echo "Unknown: $1"; show_help ;;
    esac
done

mkdir -p "$LOGDIR" "$OUTDIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/pipeline.log"; }

# ---- 准备参考序列 ----
if [ ! -f "$REF_FASTA" ]; then
    if [ -d "$VIRUS_DIR" ]; then
        log "Building reference FASTA from $VIRUS_DIR..."
        cat "$VIRUS_DIR"/*.fasta > "$REF_FASTA"
    else
        log "WARNING: $REF_FASTA not found and $VIRUS_DIR not available"
    fi
fi

# ---- 自动扫描 ----
ALL_TASKS=()
for gdir in "$SIMDIR"/group_*/; do
    [ ! -d "$gdir" ] && continue
    g=$(basename "$gdir" | sed 's/group_//')
    for mdir in "$gdir"Dataset_Mut_*pct/; do
        [ ! -d "$mdir" ] && continue
        mut=$(basename "$mdir" | sed 's/Dataset_Mut_//' | sed 's/pct$//')
        NAME="group${g}_mut${mut}"
        DONE="${OUTDIR}/${NAME}/.DONE"
        [ -f "$DONE" ] && { log "[$NAME] Skip (.DONE)"; continue; }
        INPUT="${mdir}Jackknife_Subsamples"
        [ ! -d "$INPUT" ] && { log "[$NAME] No Jackknife_Subsamples"; continue; }
        ALL_TASKS+=("$NAME|$INPUT")
    done
done

TOTAL=${#ALL_TASKS[@]}
if [ "$TOTAL" -eq 0 ]; then
    log "No pending tasks. All done."
    exit 0
fi
log "Scanned: $TOTAL pending tasks (batch=$BATCH, threads=$THREADS, jobs=$JOBS)"
log "Input:  $SIMDIR"
log "Output: $OUTDIR"
for t in "${ALL_TASKS[@]}"; do log "  - ${t%%|*}"; done

# ---- 分批运行 ----
for ((i=0; i<TOTAL; i+=BATCH)); do
    N=$((i/BATCH + 1))
    log "=== Batch $N ==="
    for ((j=i; j<i+BATCH && j<TOTAL; j++)); do
        IFS='|' read -r NAME INPUT <<< "${ALL_TASKS[$j]}"
        mkdir -p "${OUTDIR}/${NAME}"
        log "  [$NAME] Start"
        (
            t0=$(date +%s)
            /usr/bin/time -v -o "${LOGDIR}/${NAME}.time" \
                python ~/bin/assembly_pipeline.py \
                    -t all -i "$INPUT" -l 200 \
                    -o "${OUTDIR}/${NAME}/" \
                    --tmp-dir "/tmp/asm_${NAME}" \
                    --jobs "$JOBS" --threads "$THREADS" \
                    --refineC_split --refineC_merge --refineC_min_id 0.95 \
                > "${LOGDIR}/${NAME}.log" 2>&1
            ret=$?
            t1=$(date +%s)
            mem=$(grep "Maximum resident set size" "${LOGDIR}/${NAME}.time" 2>/dev/null | awk '{print $NF}')
            if [ $ret -eq 0 ]; then
                touch "${OUTDIR}/${NAME}/.DONE"
                log "  [$NAME] OK | Time: $((t1-t0))s | RSS: ${mem:-N/A}KB"
            else
                log "  [$NAME] FAIL (exit: $ret)"
            fi
        ) &
    done
    wait
    log "=== Batch $N done ==="
done

log "=== All assemblies done. Run MetaQUAST separately. ==="
