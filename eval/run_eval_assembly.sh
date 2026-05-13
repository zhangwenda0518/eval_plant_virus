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
METAQUAST_DIR="step6_metaquast_7way"
BATCH=2
THREADS=15
JOBS=4
MERGE_THREADS=60
MERGE_JOBS=5
MIN_ID=0.95
MIN_COV=0.50
FRAG_MIN_LEN=1000
SKIP_MERGE=false
SKIP_METAQUAST=false
SKIP_CHIMERIC=false

show_help() {
    echo "Usage: bash run_eval_assembly.sh [options]"
    echo ""
    echo "  阶段1: 自动扫描并分批运行 assembly_pipeline.py"
    echo "  阶段2: 补充4组 RefineC merge (MH merge / MH split+merge / ALL merge / ALL split+merge)"
    echo "  阶段3: MetaQUAST 7组对比"
    echo ""
    echo "Options:"
    echo "  --sim-data DIR      模拟数据根目录 (default: step2_benchmark_data)"
    echo "  --output-dir DIR    组装输出目录 (default: step6_assemblies)"
    echo "  --logdir DIR        日志目录 (default: step6_logs)"
    echo "  --virus-dir DIR     病毒参考基因组目录 (default: step1_eval_viruses)"
    echo "  --ref-fasta FILE    合并后的参考FASTA (default: step6_ref_viruses.fasta)"
    echo "  --metaquast-dir DIR MetaQUAST输出目录 (default: step6_metaquast_7way)"
    echo "  --batch N           每批并发组数 (default: 2)"
    echo "  --threads N         单样本线程数 (default: 15)"
    echo "  --jobs N            assembly_pipeline 内部并发 (default: 4)"
    echo "  --merge-threads N   RefineC merge 线程 (default: 60)"
    echo "  --min-id FLOAT      RefineC merge --min-id (default: 0.99)"
    echo "  --min-cov FLOAT     RefineC merge --min-cov (default: 0.90)"
    echo "  --skip-7way         跳过7组对比阶段"
    echo "  -h, --help          显示帮助"
    exit 0
}

# ---- 解析参数 ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --sim-data)    SIMDIR="$2"; shift 2 ;;
        --output-dir)  OUTDIR="$2"; shift 2 ;;
        --logdir)      LOGDIR="$2"; shift 2 ;;
        --virus-dir)   VIRUS_DIR="$2"; shift 2 ;;
        --ref-fasta)   REF_FASTA="$2"; shift 2 ;;
        --metaquast-dir) METAQUAST_DIR="$2"; shift 2 ;;
        --batch)       BATCH="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --jobs)        JOBS="$2"; shift 2 ;;
        --merge-threads) MERGE_THREADS="$2"; shift 2 ;;
        --merge-jobs)   MERGE_JOBS="$2"; shift 2 ;;
        --min-id)      MIN_ID="$2"; shift 2 ;;
        --min-cov)     MIN_COV="$2"; shift 2 ;;
        --skip-7way)   SKIP_METAQUAST=true; shift ;;
        -h|--help)     show_help ;;
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
                    --keep-temp --log_dirs "${LOGDIR}" \
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

log "=== Phase 1: All assemblies done. ==="
if [ "$SKIP_MERGE" = true ]; then
    log "Skipping RefineC merge. Run MetaQUAST manually."
    exit 0
fi

# ===== 阶段2: 补充4组 RefineC merge（分批并行）=====
log "=== Phase 2: 7-way RefineC merge (jobs=$MERGE_JOBS) ==="

# 收集所有需处理的样本
MERGE_TASKS=()
for gdir in "$OUTDIR"/group*_mut*/; do
    [ ! -d "$gdir" ] && continue
    for sdir in "$gdir"Master_*/; do
        [ ! -d "$sdir" ] && continue
        SDIR="$sdir"
        SNAME=$(basename "$sdir")
        DONE_MARK="$SDIR/.merge_7way_done"
        [ -f "$DONE_MARK" ] && continue
        MERGE_TASKS+=("$SDIR|$SNAME")
    done
done

log "Merge tasks: ${#MERGE_TASKS[@]}"
if [ ${#MERGE_TASKS[@]} -eq 0 ]; then
    log "No pending merge tasks."
else
    for ((i=0; i<${#MERGE_TASKS[@]}; i+=MERGE_JOBS)); do
        N=$((i/MERGE_JOBS + 1))
        log "=== Merge Batch $N ==="
        for ((j=i; j<i+MERGE_JOBS && j<${#MERGE_TASKS[@]}; j++)); do
            IFS='|' read -r SDIR SNAME <<< "${MERGE_TASKS[$j]}"
            DONE_MARK="$SDIR/.merge_7way_done"
            (
                MEGAHIT_CTG=$(ls "$SDIR/${SNAME}_megahit.contig.fasta"* 2>/dev/null | head -1)
                RNAVIRAL_CTG=$(ls "$SDIR/${SNAME}_rnaviralspades.contig.fasta"* 2>/dev/null | head -1)
                PENGUIN_CTG=$(ls "$SDIR/${SNAME}_penguin.contig.fasta"* 2>/dev/null | head -1)

                [ -z "$MEGAHIT_CTG" ] || [ -z "$RNAVIRAL_CTG" ] && { touch "$DONE_MARK"; exit 0; }

                # M+H merge
                MH_OUT="$SDIR/MH_merge"
                [ ! -f "$MH_OUT/merged.fasta" ] && {
                    mkdir -p "$MH_OUT"
                    refineC merge --threads 4 --contigs "$MEGAHIT_CTG" "$RNAVIRAL_CTG" \
                        --prefix MH-merge --output "$MH_OUT" --min-id "$MIN_ID" --min-cov "$MIN_COV" \
                        > "$MH_OUT/log.txt" 2>&1
                } &

                # M+H split+merge
                MH_SPLIT_OUT="$SDIR/MH_split_merge"
                [ ! -f "$MH_SPLIT_OUT/merged.fasta" ] && {
                    mkdir -p "$MH_SPLIT_OUT"
                    M_SPLIT=$(ls "$SDIR/${SNAME}_megahit_refineC"/*.split.fasta.gz 2>/dev/null | head -1)
                    H_SPLIT=$(ls "$SDIR/${SNAME}_rnaviralspades_refineC"/*.split.fasta.gz 2>/dev/null | head -1)
                    MH_SPLIT_INPUTS=""
                    [ -n "$M_SPLIT" ] && MH_SPLIT_INPUTS="$MH_SPLIT_INPUTS $M_SPLIT"
                    [ -n "$H_SPLIT" ] && MH_SPLIT_INPUTS="$MH_SPLIT_INPUTS $H_SPLIT"
                    [ -n "$MH_SPLIT_INPUTS" ] && {
                        refineC merge --threads 4 --contigs $MH_SPLIT_INPUTS \
                            --prefix MH-split-merge --output "$MH_SPLIT_OUT" \
                            --min-id "$MIN_ID" --min-cov "$MIN_COV" > "$MH_SPLIT_OUT/log.txt" 2>&1
                    }
                } &

                # ALL merge
                ALLM_OUT="$SDIR/ALL_merge"
                ALLM_CTGS="$MEGAHIT_CTG $RNAVIRAL_CTG"
                [ -n "$PENGUIN_CTG" ] && ALLM_CTGS="$ALLM_CTGS $PENGUIN_CTG"
                [ ! -f "$ALLM_OUT/merged.fasta" ] && {
                    mkdir -p "$ALLM_OUT"
                    refineC merge --threads 4 --contigs $ALLM_CTGS \
                        --prefix ALL-merge --output "$ALLM_OUT" \
                        --min-id "$MIN_ID" --min-cov "$MIN_COV" > "$ALLM_OUT/log.txt" 2>&1
                } &

                wait
                touch "$DONE_MARK"
            ) &
            log "  [$SNAME] started"
        done
        wait
        log "=== Merge Batch $N done ==="
    done
fi

log "=== Phase 2 done. Run MetaQUAST separately. ==="
