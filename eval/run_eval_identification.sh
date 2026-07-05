#!/bin/bash
# ============================================================
# 评估四：候选病毒鉴定策略比较 — 统一调度脚本
# virus_identification13.py 已内置多工具联合+对抗验证+过滤
#
# 用法: bash run_eval_identification.sh [options]
# ============================================================
set -e
BIN_DIR="$(cd "$(dirname "$0")/.." && pwd)/deps"

LOGDIR="step7_logs"
OUTDIR="step7_identification"
EVAL_SEQS="step3_identification_eval/evaluation_sequences.fasta"
LABELS="step3_identification_eval/sequence_labels.tsv"
RVDB="$HOME/database/virus-db/RVDB-v30/U-RVDBv30.0-prot.dmnd"
NR_DB="$HOME/database/nr_db/nr"
UNIREF90="$HOME/database/uniport_db/uniref90/uniref90.dmnd"
VIROIDS_DB="$HOME/database/virus-db/viroids-db/viroids.fasta.blast.db"
TAXID_FILE="$HOME/database/virus-db/taxIDs/viral_taxIDs.txt"
VIRBOT="$HOME/biosoft/virus/VirBot/VirBot.py"
THREADS=20
JOBS=5

show_help() {
    echo "Usage: bash run_eval_identification.sh [options]"
    echo ""
    echo "  运行 virus_identification13.py 对评估序列集进行多工具联合鉴定"
    echo "  (已内置 RVDB/UniRef90/VirBot/Viroids + NR对抗验证 + 多维过滤)"
    echo ""
    echo "Options:"
    echo "  --eval-seqs FILE    评估序列FASTA"
    echo "  --output-dir DIR    输出目录 (default: step7_identification)"
    echo "  --logdir DIR        日志目录 (default: step7_logs)"
    echo "  --rvdb FILE         RVDB数据库"
    echo "  --nr-db FILE        NR数据库"
    echo "  --uniref90 FILE     UniRef90数据库"
    echo "  --viroids-db FILE   类病毒BLAST数据库"
    echo "  --virbot FILE       VirBot路径"
    echo "  --threads N         线程数 (default: 20)"
    echo "  --jobs N            并发 (default: 5)"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --eval-seqs)   EVAL_SEQS="$2"; shift 2 ;;
        --output-dir)  OUTDIR="$2"; shift 2 ;;
        --logdir)      LOGDIR="$2"; shift 2 ;;
        --rvdb)        RVDB="$2"; shift 2 ;;
        --nr-db)       NR_DB="$2"; shift 2 ;;
        --uniref90)    UNIREF90="$2"; shift 2 ;;
        --viroids-db)  VIROIDS_DB="$2"; shift 2 ;;
        --virbot)      VIRBOT="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --jobs)        JOBS="$2"; shift 2 ;;
        -h|--help)     show_help ;;
        *) echo "Unknown: $1"; show_help ;;
    esac
done

mkdir -p "$LOGDIR" "$OUTDIR"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/pipeline.log"; }

DONE="${OUTDIR}/.DONE"
if [ -f "$DONE" ]; then
    log "Already completed ($DONE exists), skipping"
    exit 0
fi

log "Running virus_identification13.py on $EVAL_SEQS"
t0=$(date +%s)
/usr/bin/time -v -o "$LOGDIR/identify.time" \
    python "$BIN_DIR/virus_identification.py" \
        --input "$EVAL_SEQS" -ext .fasta \
        --output "$OUTDIR/" \
        --identify_tools all --blast_mode filter \
        --virbot_path "$VIRBOT" \
        --virus_taxid "$TAXID_FILE" \
        --virus_protein_db "$RVDB" \
        --uniprot_db "$UNIREF90" \
        --viroids_db "$VIROIDS_DB" \
        --nr_db "$NR_DB" \
        --jobs "$JOBS" --threads "$THREADS" \
    > "$LOGDIR/identify.log" 2>&1
ret=$?
t1=$(date +%s)
mem=$(grep "Maximum resident set size" "$LOGDIR/identify.time" 2>/dev/null | awk '{print $NF}')

if [ $ret -eq 0 ]; then
    touch "$DONE"
    log "Done | Time: $((t1-t0))s | RSS: ${mem:-N/A}KB"
else
    log "FAIL (exit: $ret) | Time: $((t1-t0))s"
    exit $ret
fi
