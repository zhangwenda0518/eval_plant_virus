#!/bin/bash
# ============================================================
# 评估四：病毒分类方法比较 — 统一调度脚本 v2
# 合并所有测试序列为单文件，一次运行（而非500次）
#
# 用法: bash run_eval_classification.sh [options]
# ============================================================
set -e

LOGDIR="step9_logs"
OUTDIR="step9_classification"
TEST_DIR="step4_classification_eval/test_sequences"
META="step4_classification_eval/test_metadata.tsv"
MMSEQS_DB="$HOME/database/virus-db/RVDB-v31/RVDB.mmseqs_db/RVDB.mmseqs"
VITAP_DB="$HOME/database/virus-db/vitap-db/VMR-MSL40_DB"
ACVIRUS_DB="$HOME/database/virus-db/acvirus_db"
THREADS=32
JOBS=3

show_help() {
    echo "Usage: bash run_eval_classification.sh [options]"
    echo ""
    echo "  合并 test_sequences/ 下所有 .fasta 为单文件，一次运行三种分类器 + 整合"
    echo ""
    echo "Options:"
    echo "  --test-dir DIR      测试序列目录 (default: step4_classification_eval/test_sequences)"
    echo "  --meta FILE         测试真值 TSV"
    echo "  --output-dir DIR    输出目录"
    echo "  --mmseqs-db DIR     MMseqs2 数据库"
    echo "  --vitap-db DIR      VITAP 数据库"
    echo "  --acvirus-db DIR    ACVirus 数据库"
    echo "  --threads N         线程 (default: 32)"
    echo "  --jobs N            并发 (default: 3)"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --test-dir)   TEST_DIR="$2"; shift 2 ;;
        --meta)       META="$2"; shift 2 ;;
        --output-dir) OUTDIR="$2"; shift 2 ;;
        --logdir)     LOGDIR="$2"; shift 2 ;;
        --mmseqs-db)  MMSEQS_DB="$2"; shift 2 ;;
        --vitap-db)   VITAP_DB="$2"; shift 2 ;;
        --acvirus-db) ACVIRUS_DB="$2"; shift 2 ;;
        --threads)    THREADS="$2"; shift 2 ;;
        --jobs)       JOBS="$2"; shift 2 ;;
        -h|--help)    show_help ;;
        *) echo "Unknown: $1"; show_help ;;
    esac
done

mkdir -p "$LOGDIR" "$OUTDIR"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/pipeline.log"; }

# ---- 合并测试序列到单个目录 ----
MERGED_DIR="${OUTDIR}/merged_seqs"
MERGED="${MERGED_DIR}/test_all.fasta"
if [ ! -f "$MERGED" ]; then
    log "Merging test sequences from $TEST_DIR into $MERGED_DIR ..."
    mkdir -p "$MERGED_DIR"
    cat "$TEST_DIR"/*.fasta > "$MERGED"
    N=$(grep -c "^>" "$MERGED")
    log "Merged $N sequences"
else
    N=$(grep -c "^>" "$MERGED")
    log "Merged exists: $N sequences"
fi

# ---- 三工具分类 ----
CLASSIFY_DIR="${OUTDIR}/classify"
DONE="${CLASSIFY_DIR}/.DONE"
if [ -f "$DONE" ]; then
    log "[classify] Skip (.DONE)"
else
    log "[classify] Running MMseqs2 + VITAP + ACVirus on merged FASTA..."
    t0=$(date +%s)
    /usr/bin/time -v -o "$LOGDIR/classify.time" \
        python ~/bin/virus_classifier.py \
            -i "$MERGED_DIR" --ext .fasta \
            -j "$JOBS" -p "$THREADS" \
            -t mmseqs,acvirus,vitap \
            --output-dir "$CLASSIFY_DIR/" \
            --vitap-db "$VITAP_DB" \
            --mmseqs-db "$MMSEQS_DB" \
            --acvirus-db "$ACVIRUS_DB" \
        > "$LOGDIR/classify.log" 2>&1
    ret=$?; t1=$(date +%s)
    mem=$(grep "Maximum resident set size" "$LOGDIR/classify.time" 2>/dev/null | awk '{print $NF}')
    if [ $ret -eq 0 ]; then
        mkdir -p "$CLASSIFY_DIR"; touch "$DONE"
        log "[classify] OK | Time: $((t1-t0))s | RSS: ${mem:-N/A}KB"
    else
        log "[classify] FAIL (exit: $ret)"
        exit $ret
    fi
fi

# ---- 整合分类结果 ----
INTEGRATE_DIR="${OUTDIR}/integrated"
DONE="${INTEGRATE_DIR}/.DONE"
if [ -f "$DONE" ]; then
    log "[integrate] Skip (.DONE)"
else
    MM=$(ls "${CLASSIFY_DIR}/mmseqs_results/"*_lca.tsv 2>/dev/null | head -1)
    VT=$(ls "${CLASSIFY_DIR}/VITAP_results/"*.best_determined_lineages.tsv 2>/dev/null | head -1)
    AC=$(ls "${CLASSIFY_DIR}/ACVirus_results/"*.final_result.tsv 2>/dev/null | head -1)

    if [ -z "$MM" ] || [ -z "$VT" ] || [ -z "$AC" ]; then
        log "[integrate] Missing results (MM=$MM VT=$VT AC=$AC)"
        exit 1
    fi

    log "[integrate] Running virus_classifier_analysis13.R ..."
    t0=$(date +%s)
    Rscript ~/bin/virus_classifier_analysis13.R \
        --mmseqs "$MM" \
        --vitap "$VT" \
        --acvirus "$AC" \
        -o "$INTEGRATE_DIR/" \
        > "$LOGDIR/integrate.log" 2>&1
    ret=$?; t1=$(date +%s)
    if [ $ret -eq 0 ]; then
        mkdir -p "$INTEGRATE_DIR"; touch "$DONE"
        log "[integrate] OK | Time: $((t1-t0))s"
    else
        log "[integrate] FAIL (exit: $ret)"
        exit $ret
    fi
fi

log "Done."
