#!/bin/bash
# ============================================================
# 评估二：已知病毒检测方法比较 — 统一调度脚本 v2
# 支持断点续跑、时间/内存记录、预构建索引自动探测
# ============================================================
set -e
BIN_DIR="$(cd "$(dirname "$0")/.." && pwd)/deps"

# ---- 默认值 ----
LOGDIR="step5_logs"
REF_INFO="final.cluster.ref_info.tsv"
REF_FASTA="final.cluster.ref.fasta"
REF_DB=""            # 可选：预构建数据库根目录
GENES_COV="$HOME/database/virus-db/plant_virus_db/ref_db/virus_genes_cov.tsv"
INPUT_DIR="step5_LoD_all"
OUTDIR="step5_results"
MULTITOOL_DB="$HOME/database/virus-db/plantvirus-db/plant_virus_db/5.virus.ref.build.db"
THREADS=10
ALIGN_THREADS=8
MULTITOOL_THREADS=20
MULTITOOL_JOBS=2
SKIP_SALMON=false
SKIP_BOWTIE2=false
SKIP_MINIMAP2=false
SKIP_MULTITOOL=false

# ---- 解析命名参数 ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --logdir)          LOGDIR="$2"; shift 2 ;;
        --ref-info)        REF_INFO="$2"; shift 2 ;;
        --ref-fasta)       REF_FASTA="$2"; shift 2 ;;
        --ref-db)          REF_DB="$2"; shift 2 ;;
        --genes-cov)       GENES_COV="$2"; shift 2 ;;
        --input-dir)       INPUT_DIR="$2"; shift 2 ;;
        --output-dir)      OUTDIR="$2"; shift 2 ;;
        --multitool-db)    MULTITOOL_DB="$2"; shift 2 ;;
        --threads)         THREADS="$2"; shift 2 ;;
        --align-threads)   ALIGN_THREADS="$2"; shift 2 ;;
        --multitool-threads) MULTITOOL_THREADS="$2"; shift 2 ;;
        --multitool-jobs)  MULTITOOL_JOBS="$2"; shift 2 ;;
        --tools)           MULTITOOL_TOOLS="$2"; shift 2 ;;
        --skip-salmon)     SKIP_SALMON=true; shift ;;
        --skip-bowtie2)    SKIP_BOWTIE2=true; shift ;;
        --skip-minimap2)   SKIP_MINIMAP2=true; shift ;;
        --skip-multitool)  SKIP_MULTITOOL=true; shift ;;
        -h|--help)
            echo "Usage: bash run_eval_known_virus.sh [options]"
            echo ""
            echo "Options:"
            echo "  --logdir DIR            日志目录 (default: step5_logs)"
            echo "  --ref-fasta FILE        参考FASTA (default: final.cluster.ref.fasta)"
            echo "  --ref-info FILE         ref_info TSV"
            echo "  --ref-db DIR            预构建数据库根目录（含 index_db/kraken2_db/等）"
            echo "  --genes-cov FILE        基因覆盖TSV"
            echo "  --input-dir DIR         LoD数据目录 (default: step5_LoD_all)"
            echo "  --output-dir DIR        输出根目录 (default: step5_results)"
            echo "  --multitool-db DIR      预构建数据库根目录 (含kraken2_db/kraken2x_db等)"
            echo "  --threads N             并发样本数 (default: 10)"
            echo "  --align-threads N       单样本比对线程 (default: 8)"
            echo "  --multitool-threads N   batch_class 线程 (default: 20)"
            echo "  --tools LIST            batch_class 工具列表 (逗号分隔)"
            echo "  --skip-salmon           跳过 Salmon"
            echo "  --skip-bowtie2          跳过 Bowtie2"
            echo "  --skip-minimap2         跳过 Minimap2"
            echo "  --skip-multitool        跳过批量多工具"
            exit 0
            ;;
        *) echo "Unknown: $1"; exit 1 ;;
    esac
done

MULTITOOL_TOOLS="${MULTITOOL_TOOLS:-kraken2,kraken2x,krakenuniq,centrifuger,kunpeng,metabuli,sylph,kaiju}"

mkdir -p "$LOGDIR" "$OUTDIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/pipeline.log"
}

# ---- 合并 LoD 数据 ----
if [ ! -d "$INPUT_DIR" ] || [ "$(ls "$INPUT_DIR"/*_R1.fastq.gz 2>/dev/null | wc -l)" -eq 0 ]; then
    log "Merging LoD data..."
    mkdir -p "$INPUT_DIR"
    for g in 1 2 3 4 5; do
        SRC="step2_benchmark_data/group_$g/Dataset_LoD_Test"
        if [ -d "$SRC" ]; then
            for f in "$SRC"/LoD_Mixed_*_R1.fastq.gz; do
                [ -e "$f" ] && ln -sf "$(realpath "$f")" "$INPUT_DIR/" 2>/dev/null || true
            done
            for f in "$SRC"/LoD_Mixed_*_R2.fastq.gz; do
                [ -e "$f" ] && ln -sf "$(realpath "$f")" "$INPUT_DIR/" 2>/dev/null || true
            done
        fi
    done
    log "Total LoD samples: $(ls "$INPUT_DIR"/*_R1.fastq.gz 2>/dev/null | wc -l)"
fi

# ---- 索引预处理：从预构建数据库探测并链接已有索引 ----
setup_index() {
    local tool="$1"
    local outdir="$2"

    # 如果没有指定预构建数据库，跳过（让 batch_virus_depth40 自己构建）
    if [ -z "$REF_DB" ] || [ ! -d "$REF_DB" ]; then
        return 0
    fi

    local idx_dir="$outdir/index"
    mkdir -p "$idx_dir"

    local fasta_name
    fasta_name=$(basename "$REF_FASTA")

    case "$tool" in
        salmon)
            local src="$REF_DB/ref.virus.build.index_db/ref.virus.build.sal_index"
            # 也检查可能的其他命名
            [ ! -e "$src" ] && src="$REF_DB/ref.virus.build.index_db/ref.virus.build.sal_index"
            local dst="$idx_dir/${fasta_name}.salmon_idx"
            ;;
        bowtie2)
            # Bowtie2 索引是一组 .bt2 文件，按前缀查找
            local src_prefix="$REF_DB/ref.virus.build.index_db/ref.virus.build.bt2_index"
            # 检查是否有提前建好的 bowtie2 索引
            if ls "$REF_DB"/ref.virus.build.index_db/*.1.bt2 2>/dev/null | head -1 > /dev/null; then
                local bt2_base=$(ls "$REF_DB"/ref.virus.build.index_db/*.1.bt2 2>/dev/null | head -1 | sed 's/\.1\.bt2//')
                for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
                    if [ -f "${bt2_base}.${ext}" ]; then
                        ln -sf "${bt2_base}.${ext}" "$idx_dir/${fasta_name}_bowtie2.${ext}"
                    fi
                done
            fi
            return 0  # Bowtie2 索引结构复杂，探测成功就链接，否则让脚本自建
            ;;
        minimap2)
            # Minimap2 用 FASTA 本身，不需要索引
            return 0
            ;;
        *)
            return 0
            ;;
    esac

    if [ -e "$src" ] && [ ! -e "$dst" ]; then
        ln -sf "$src" "$dst"
        log "[$tool] Linked prebuilt index: $src -> $dst"
    elif [ -e "$dst" ]; then
        log "[$tool] Index already exists: $dst"
    fi
}

# ---- 运行函数 ----
run_with_log() {
    local name="$1"; shift
    local outdir="$1"; shift
    local logfile="$LOGDIR/${name}.log"

    if [ -f "$outdir/.DONE" ]; then
        log "[$name] Already completed, skipping"
        return 0
    fi

    # 预建索引链接
    setup_index "$name" "$outdir"

    log "[$name] Starting... (log: $logfile)"
    local t0
    t0=$(date +%s)
    /usr/bin/time -v -o "$LOGDIR/${name}.time" "$@" > "$logfile" 2>&1
    local ret=$?
    local t1
    t1=$(date +%s)
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

# ---- 1. Salmon ----
if [ "$SKIP_SALMON" != true ]; then
    run_with_log "salmon" "${OUTDIR}/salmon" \
        python "$BIN_DIR/batch_virus_depth.py" \
            --input_dir "$INPUT_DIR" \
            --output_dir "${OUTDIR}/salmon" \
            --ref_info "$REF_INFO" \
            --reference "$REF_FASTA" \
            --genes_cov "$GENES_COV" \
            --tool salmon --threads "$THREADS" --align_threads "$ALIGN_THREADS" \
            --batch_size 100 --resume || true
fi

# ---- 2. Bowtie2 ----
if [ "$SKIP_BOWTIE2" != true ]; then
    run_with_log "bowtie2" "${OUTDIR}/bowtie2" \
        python "$BIN_DIR/batch_virus_depth.py" \
            --input_dir "$INPUT_DIR" \
            --output_dir "${OUTDIR}/bowtie2" \
            --ref_info "$REF_INFO" \
            --reference "$REF_FASTA" \
            --genes_cov "$GENES_COV" \
            --tool bowtie2 --threads "$THREADS" --align_threads "$ALIGN_THREADS" \
            --batch_size 100 --resume || true
fi

# ---- 3. Minimap2 ----
if [ "$SKIP_MINIMAP2" != true ]; then
    run_with_log "minimap2" "${OUTDIR}/minimap2" \
        python "$BIN_DIR/batch_virus_depth.py" \
            --input_dir "$INPUT_DIR" \
            --output_dir "${OUTDIR}/minimap2" \
            --ref_info "$REF_INFO" \
            --reference "$REF_FASTA" \
            --genes_cov "$GENES_COV" \
            --tool minimap2 --threads "$THREADS" --align_threads "$ALIGN_THREADS" \
            --batch_size 100 --resume || true
fi

# ---- 4. 批量多工具 ----
if [ "$SKIP_MULTITOOL" != true ]; then
    run_with_log "multitool" "${OUTDIR}/multitool" \
        python "$BIN_DIR/batch_class.reads.py" \
            -i "$INPUT_DIR" \
            -o "${OUTDIR}/multitool" \
            --db-dir "$MULTITOOL_DB" \
            --tools kraken2 kraken2x krakenuniq centrifuger kunpeng metabuli sylph kaiju \
            --jobs "$MULTITOOL_JOBS" --threads "$MULTITOOL_THREADS" || true
fi

# ---- 汇总 ----
log "====================="
log "Evaluation Summary"
log "====================="
for name in salmon bowtie2 minimap2 multitool; do
    if [ -f "$LOGDIR/${name}.time" ]; then
        elapsed=$(grep "Elapsed" "$LOGDIR/${name}.time" 2>/dev/null | awk '{print $NF}' || echo "N/A")
        rss=$(grep "Maximum resident set size" "$LOGDIR/${name}.time" 2>/dev/null | awk '{print $NF}' || echo "N/A")
        if [ -f "${OUTDIR}/${name}/.DONE" ]; then
            echo "✅ $name | Time=$elapsed | Peak RSS=${rss}KB"
        else
            echo "❌ $name | Time=$elapsed | Peak RSS=${rss}KB"
        fi
    fi
done | tee -a "$LOGDIR/pipeline.log"

log "Pipeline finished."
