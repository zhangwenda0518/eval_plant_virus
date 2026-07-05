#!/bin/bash
# ============================================================
# 评估二补充: MetaQUAST 7组 + 嵌合体检测
# 用法: bash run_assembly_metaquast.sh [ASMDIR] [REF] [OUTDIR]
# ============================================================
set -e
BIN_DIR="$(cd "$(dirname "$0")/.." && pwd)/deps"

ASMDIR="${1:-step6_assemblies}"
REF="${2:-step6_ref_viruses.fasta}"
OUTDIR="${3:-step6_metaquast_7way}"
THREADS="${4:-4}"

mkdir -p "$OUTDIR"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

for gdir in "$ASMDIR"/group*_mut*/; do
    [ ! -d "$gdir" ] && continue
    for sdir in "$gdir"1.virus-assembly/*/; do
        [ ! -d "$sdir" ] && continue
        SDIR="$sdir"
        SNAME=$(basename "$sdir")
        QUAST_OUT="$OUTDIR/$SNAME"
        DONE="$QUAST_OUT/.DONE"
        [ -f "$DONE" ] && { log "[$SNAME] skip"; continue; }
        mkdir -p "$QUAST_OUT"

        CTGS=""; LABELS=""
        MF=$(ls "$SDIR/${SNAME}_megahit.contig.fasta"* 2>/dev/null | head -1)
        RF=$(ls "$SDIR/${SNAME}_rnaviralspades.contig.fasta"* 2>/dev/null | head -1)
        PF=$(ls "$SDIR/${SNAME}_penguin.contig.fasta"* 2>/dev/null | head -1)
        MHF="$SDIR/MH_merge/merged.fasta"
        MHSF="$SDIR/MH_split_merge/merged.fasta"
        AMF="$SDIR/ALL_merge/merged.fasta"
        ASMF=$(ls "$SDIR/${SNAME}_all_tools_refineC_merge.merged.fasta"* 2>/dev/null | head -1)

        [ -n "$MF" ]   && CTGS="$CTGS $MF"   && LABELS="$LABELS Megahit"
        [ -n "$RF" ]   && CTGS="$CTGS $RF"   && LABELS="$LABELS RNAViralSPAdes"
        [ -n "$PF" ]   && CTGS="$CTGS $PF"   && LABELS="$LABELS Penguin"
        [ -f "$MHF" ]  && CTGS="$CTGS $MHF"  && LABELS="$LABELS MH_Merge"
        [ -f "$MHSF" ] && CTGS="$CTGS $MHSF" && LABELS="$LABELS MH_SplitMerge"
        [ -f "$AMF" ]  && CTGS="$CTGS $AMF"  && LABELS="$LABELS ALL_Merge"
        [ -f "$ASMF" ] && CTGS="$CTGS $ASMF" && LABELS="$LABELS ALL_SplitMerge"

        [ -z "$CTGS" ] && { log "[$SNAME] No contigs"; continue; }
        metaquast -o "$QUAST_OUT" -r "$REF" \
            --min-contig 200 -t "$THREADS" -l "$LABELS" $CTGS \
            > "$QUAST_OUT/metaquast.log" 2>&1 &
        touch "$DONE"
    done
done
wait
log "MetaQUAST done: $OUTDIR"

# ===== 嵌合体检测 =====
log "=== Chimeric contig detection ==="
CHIMDIR="${OUTDIR}_chimeric"
mkdir -p "$CHIMDIR"

for gdir in "$ASMDIR"/group*_mut*/; do
    [ ! -d "$gdir" ] && continue
    for sdir in "$gdir"1.virus-assembly/*/; do
        SDIR="$sdir"
        SNAME=$(basename "$sdir")

        for ctg_label in megahit rnaviralspades penguin MH_merge MH_split_merge ALL_merge all_tools_refineC_merge; do
            case $ctg_label in
                megahit) CTG=$(ls "$SDIR/${SNAME}_megahit.contig.fasta"* 2>/dev/null | head -1) ;;
                rnaviralspades) CTG=$(ls "$SDIR/${SNAME}_rnaviralspades.contig.fasta"* 2>/dev/null | head -1) ;;
                penguin) CTG=$(ls "$SDIR/${SNAME}_penguin.contig.fasta"* 2>/dev/null | head -1) ;;
                MH_merge) CTG="$SDIR/MH_merge/merged.fasta" ;;
                MH_split_merge) CTG="$SDIR/MH_split_merge/merged.fasta" ;;
                ALL_merge) CTG="$SDIR/ALL_merge/merged.fasta" ;;
                all_tools_refineC_merge) CTG=$(ls "$SDIR/${SNAME}_all_tools_refineC_merge.merged.fasta"* 2>/dev/null | head -1) ;;
            esac
            [ ! -f "$CTG" ] && continue
            OUT="$CHIMDIR/${SNAME}_${ctg_label}_blastn.tsv"
            RPT="$CHIMDIR/${SNAME}_${ctg_label}_chimeric.tsv"
            [ -f "$RPT" ] && continue
            blastn -query "$CTG" -db "$REF" \
                -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
                -out "$OUT" -evalue 1e-5 -num_threads 4 &
        done
    done
done
wait

# 汇总嵌合率
for blast_out in "$CHIMDIR"/*_blastn.tsv; do
    [ ! -f "$blast_out" ] && continue
    python "$BIN_DIR/detect_chimeric_contigs.py" \
        -i "$blast_out" -o "${blast_out%.tsv}_chimeric.tsv" 2>/dev/null
done
log "Chimeric detection done: $CHIMDIR"
