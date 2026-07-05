#!/bin/bash
# ============================================================
# 评估七：宿主分类基准 — 三工具集成 + 决策树 + 评估
#
# RNAVirHost + PhaBOX2 + ICTV(C9) → 决策树集成 → 评估
#
# 用法:
#   bash run_eval_host_prediction.sh
#   bash run_eval_host_prediction.sh --input contigs.fasta --tax taxonomy.tsv --mode all
# ============================================================
set -e
BIN_DIR="$(cd "$(dirname "$0")/.." && pwd)/deps"

# ---- 默认值 ----
INPUT="step9_classification/classifier/evaluation_sequences.fasta"
TAX="step9_classification/integrated/final_integrated_classification.tsv"
META="step4_classification_eval/test_metadata_full.tsv"
LABELS="step3_master_eval/sequence_labels_category.tsv"
OUTDIR="step11_host_evaluation"
PHABOX_DB="$HOME/database/virus-db/phabox_db_v2_2"
PROB_DIR="$HOME/database/virus-db/plantvirus-db/plant_virus_db/1.virus-host_db/C-host_classify/cross_analysis"
THREADS=40

# ---- 解析参数 ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --input) INPUT="$2"; shift 2 ;;
        --tax) TAX="$2"; shift 2 ;;
        --meta) META="$2"; shift 2 ;;
        --labels) LABELS="$2"; shift 2 ;;
        --output-dir) OUTDIR="$2"; shift 2 ;;
        --phabox-db) PHABOX_DB="$2"; shift 2 ;;
        --prob-dir) PROB_DIR="$2"; shift 2 ;;
        --threads|-t) THREADS="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --force|-f) FORCE="--force"; shift ;;
        --skip-rnavirhost) SKIP_RVH="--skip-rnavirhost"; shift ;;
        --skip-phabox) SKIP_PB="--skip-phabox"; shift ;;
        --skip-ictv) SKIP_ICTV="--skip-ictv"; shift ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

mkdir -p "$OUTDIR" "${OUTDIR}/evaluation"

# ---- 阶段1: 运行宿主预测 (三工具 + 决策树集成) ----
DONE="${OUTDIR}/.DONE"
if [ -f "$DONE" ] && [ -z "$FORCE" ]; then
    echo "[host_prediction] Skip (.DONE exists)"
else
    echo "[host_prediction] Running: RNAVirHost + PhaBOX2 + ICTV(C9) → Ensemble"
    /usr/bin/time -v -o "${OUTDIR}/host_prediction.time" \
        python "$BIN_DIR/run_host_prediction.py" \
            --input "$INPUT" \
            --tax "$TAX" \
            --output-dir "$OUTDIR" \
            --threads "$THREADS" \
            --phabox-db "$PHABOX_DB" \
            --prob-dir "$PROB_DIR" \
            ${MODE:+--mode "$MODE"} \
            ${SKIP_RVH} ${SKIP_PB} ${SKIP_ICTV} ${FORCE}
    touch "$DONE"
fi

# ---- 阶段2: 评估 ----
echo "[evaluation] Running eval_host_prediction.py ..."
python "$(dirname "$0")/../metrics/eval_host_prediction.py" \
    --rvh "${OUTDIR}/RVH_result/result.csv" \
    --phabox "${OUTDIR}/phabox2_output/final_prediction/cherry_prediction.tsv" \
    --c9 "${OUTDIR}/C9_ICTV_result/classification_result.tsv" \
    --ensemble "${OUTDIR}/ensemble_host_summary.tsv" \
    --labels "$LABELS" \
    --outdir "${OUTDIR}/evaluation/"

echo ""
echo "Done: $OUTDIR"
echo "  Results: ${OUTDIR}/evaluation/host_prediction_metrics.tsv"
echo "  Figures: ${OUTDIR}/evaluation/Fig_Host_*.png"
