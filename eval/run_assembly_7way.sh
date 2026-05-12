#!/bin/bash
# ============================================================
# 组装评估 7组对比 — RefineC merge 组合 + MetaQUAST
# ============================================================
set -e

SIMDIR="${1:-step2_Ultimate_Data}"
ASMDIR="${2:-step6_assemblies}"
OUTDIR="${3:-step6_metaquast_7way}"
REF="${4:-step6_ref_viruses.fasta}"
THREADS="${5:-60}"
MIN_ID="${6:-0.99}"
MIN_COV="${7:-0.90}"

mkdir -p "$OUTDIR"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# 自动扫描所有组装产出的样本目录
SAMPLES=()
for gdir in "$ASMDIR"/group*_mut*/; do
    [ ! -d "$gdir" ] && continue
    for sdir in "$gdir"1.virus-assembly/*/; do
        [ ! -d "$sdir" ] && continue
        sname=$(basename "$sdir")
        SAMPLES+=("$sdir|$sname")
    done
done

log "Found ${#SAMPLES[@]} assembly sample directories"
if [ ${#SAMPLES[@]} -eq 0 ]; then
    log "ERROR: No samples found"
    exit 1
fi

for entry in "${SAMPLES[@]}"; do
    IFS='|' read -r SDIR SNAME <<< "$entry"

    MEGAHIT_CTG="$SDIR/${SNAME}_megahit.contig.fasta"
    RNAVIRAL_CTG="$SDIR/${SNAME}_rnaviralspades.contig.fasta"
    PENGUIN_CTG="$SDIR/${SNAME}_penguin.contig.fasta"
    ALL_MERGED="$SDIR/${SNAME}_all_tools_refineC_merge.merged.fasta"

    # 找 fasta 或 fasta.gz
    for f in "$MEGAHIT_CTG" "$MEGAHIT_CTG.gz"; do [ -f "$f" ] && MEGAHIT_CTG="$f" && break; done
    for f in "$RNAVIRAL_CTG" "$RNAVIRAL_CTG.gz"; do [ -f "$f" ] && RNAVIRAL_CTG="$f" && break; done
    for f in "$PENGUIN_CTG" "$PENGUIN_CTG.gz"; do [ -f "$f" ] && PENGUIN_CTG="$f" && break; done
    for f in "$ALL_MERGED" "$ALL_MERGED.gz"; do [ -f "$f" ] && ALL_MERGED="$f" && break; done

    DONE_MARK="$SDIR/.merge_7way_done"
    [ -f "$DONE_MARK" ] && { log "[$SNAME] 7-way merge already done, skip"; continue; }

    log "[$SNAME] Running 4 RefineC merge combinations..."

    # === M+H merge (不split) ===
    MH_OUT="$SDIR/MH_merge"
    if [ ! -f "$MH_OUT/merged.fasta" ]; then
        mkdir -p "$MH_OUT"
        /usr/bin/time -v -o "$MH_OUT/time.log" \
            refineC merge --threads "$THREADS" \
                --contigs "$MEGAHIT_CTG" "$RNAVIRAL_CTG" \
                --prefix MH-merge --output "$MH_OUT" \
                --min-id "$MIN_ID" --min-cov "$MIN_COV" \
            > "$MH_OUT/log.txt" 2>&1 || log "[$SNAME] MH merge failed"
    fi

    # === M+H split+merge ===
    MH_SPLIT_OUT="$SDIR/MH_split_merge"
    if [ ! -f "$MH_SPLIT_OUT/merged.fasta" ]; then
        # 先 split
        MH_SPLIT_TMP="$SDIR/MH_split_tmp"
        mkdir -p "$MH_SPLIT_TMP" "$MH_SPLIT_OUT"
        cat "$MEGAHIT_CTG" "$RNAVIRAL_CTG" > "$MH_SPLIT_TMP/combined.fasta"
        /usr/bin/time -v -o "$MH_SPLIT_OUT/time_split.log" \
            refineC split --threads "$THREADS" \
                --contigs "$MH_SPLIT_TMP/combined.fasta" \
                --prefix MH-split --output "$MH_SPLIT_TMP/split" \
                --frag-min-len 300 \
            > "$MH_SPLIT_OUT/log_split.txt" 2>&1 || log "[$SNAME] MH split failed"
        # 再 merge
        MH_SPLIT_FA=$(ls "$MH_SPLIT_TMP/split"/*.fasta 2>/dev/null | head -1)
        [ -z "$MH_SPLIT_FA" ] && MH_SPLIT_FA="$MH_SPLIT_TMP/combined.fasta"
        /usr/bin/time -v -o "$MH_SPLIT_OUT/time_merge.log" \
            refineC merge --threads "$THREADS" \
                --contigs "$MH_SPLIT_FA" \
                --prefix MH-split-merge --output "$MH_SPLIT_OUT" \
                --min-id "$MIN_ID" --min-cov "$MIN_COV" \
            > "$MH_SPLIT_OUT/log_merge.txt" 2>&1 || log "[$SNAME] MH split+merge failed"
    fi

    # === 全三者 merge (不split) ===
    ALLM_OUT="$SDIR/ALL_merge"
    if [ ! -f "$ALLM_OUT/merged.fasta" ]; then
        mkdir -p "$ALLM_OUT"
        /usr/bin/time -v -o "$ALLM_OUT/time.log" \
            refineC merge --threads "$THREADS" \
                --contigs "$MEGAHIT_CTG" "$RNAVIRAL_CTG" "$PENGUIN_CTG" \
                --prefix ALL-merge --output "$ALLM_OUT" \
                --min-id "$MIN_ID" --min-cov "$MIN_COV" \
            > "$ALLM_OUT/log.txt" 2>&1 || log "[$SNAME] ALL merge failed"
    fi

    touch "$DONE_MARK"
    log "[$SNAME] 7-way merge done"
done

# ===== MetaQUAST 7组 =====
log "Running MetaQUAST for 7-way comparison..."

for entry in "${SAMPLES[@]}"; do
    IFS='|' read -r SDIR SNAME <<< "$entry"

    QUAST_OUT="$OUTDIR/$SNAME"
    mkdir -p "$QUAST_OUT"

    CTG_LIST=""
    LABELS=""

    # 1. MEGAHIT
    MF=$(ls "$SDIR/${SNAME}_megahit.contig.fasta"* 2>/dev/null | head -1)
    [ -n "$MF" ] && CTG_LIST="$CTG_LIST $MF" && LABELS="$LABELS Megahit"

    # 2. rnaviralSPAdes
    RF=$(ls "$SDIR/${SNAME}_rnaviralspades.contig.fasta"* 2>/dev/null | head -1)
    [ -n "$RF" ] && CTG_LIST="$CTG_LIST $RF" && LABELS="$LABELS RNAViralSPAdes"

    # 3. Penguin
    PF=$(ls "$SDIR/${SNAME}_penguin.contig.fasta"* 2>/dev/null | head -1)
    [ -n "$PF" ] && CTG_LIST="$CTG_LIST $PF" && LABELS="$LABELS Penguin"

    # 4. M+H merge
    MHF="$SDIR/MH_merge/merged.fasta"
    [ -f "$MHF" ] && CTG_LIST="$CTG_LIST $MHF" && LABELS="$LABELS MH_Merge"

    # 5. M+H split+merge
    MHSF="$SDIR/MH_split_merge/merged.fasta"
    [ -f "$MHSF" ] && CTG_LIST="$CTG_LIST $MHSF" && LABELS="$LABELS MH_SplitMerge"

    # 6. 全三者 merge
    AMF="$SDIR/ALL_merge/merged.fasta"
    [ -f "$AMF" ] && CTG_LIST="$CTG_LIST $AMF" && LABELS="$LABELS ALL_Merge"

    # 7. 全三者 split+merge (assembly_pipeline 默认产物)
    ASMF=$(ls "$SDIR/${SNAME}_all_tools_refineC_merge.merged.fasta"* 2>/dev/null | head -1)
    [ -f "$ASMF" ] && CTG_LIST="$CTG_LIST $ASMF" && LABELS="$LABELS ALL_SplitMerge"

    if [ -n "$CTG_LIST" ]; then
        metaquast -o "$QUAST_OUT" -r "$REF" \
            --min-contig 200 -t 4 \
            -l "$LABELS" $CTG_LIST \
            > "$QUAST_OUT/metaquast.log" 2>&1 &
    fi
done
wait

log "7-way MetaQUAST done. Results: $OUTDIR"
