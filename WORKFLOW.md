# 第二章：病毒相关基准分析

## 2.2.1 模拟数据集

```bash
# ====== 前置数据：参考病毒信息 ======
# final.cluster.ref_info.tsv  final.cluster.ref.fasta  (软链接)

# ====== Step 1: 选取病毒（50非节段 + 10节段） ======
python prep/prep_select_eval_viruses.py \
    --ref-info final.cluster.ref_info.tsv \
    --ref-fasta final.cluster.ref.fasta \
    --n-viruses 50 --n-segmented 10 \
    --include "NC_002030.1,OR489165.1" \
    --outdir step1_eval_viruses/ --seed 42
```

## 2.2.2 去宿主评估

### 2.2.2.1 数据模拟
```bash
# ====== Step 2: 评估1 宿主过滤 -- 60病毒混合 + 宿主背景 ======
mkdir step2_simulator
python sim/virome_simulator.py LoD_mix \
    -i step1_eval_viruses/ --bgref ningxia.genome.fasta \
    --bg-reads 10000000 --depths 1000 --mode PE -l 150 \
    --seed 42 --threads 30 -o step2_simulator/eval1_host_dep --all-in-one
```

### 2.2.2.2 宿主去除
```bash
mkdir -p step5_sim && cd step5_sim
ln -s ../step2_simulator/eval1_host_dep/Host_Depletion_Mixed_PE_R* ./
cd ..

bash eval/run_eval_host_depletion.sh \
    --host-k2 ~/database/host_db/kraken2/ \
    --host-hisat2 ~/database/host_db/hisat2/host \
    --jobs 2 --threads 20 \
    --input-dir step5_sim \
    --output-dir step5_host_free
```

### 2.2.2.3 结果分析
```bash
# 统计
python metrics/eval_host_depletion.py \
    -i step5_host_free --outdir step5_host_free_analysis/

# 比较绘图
python plots/plot_host_depletion.py \
    --input step5_host_free_analysis/host_depletion_detail.tsv \
    --out_dir step5_host_free_analysis --dpi 600

# 资源消耗绘图
python plots/plot_host_depletion.resources.py \
    --log_dir step5_host_free_logs/ \
    --out_dir step5_host_free_analysis
```

## 2.2.3 已知病毒鉴定基准分析

### 2.2.3.1 数据模拟
```bash
# ====== Step 3: 评估2A 已知检测(丰度定量) -- 5M 病毒对数正态分配 ======
python sim/virome_simulator.py gen-config \
    -i step1_eval_viruses/ -o step2_simulator/eval2a_config.txt \
    --total-reads 5M

python sim/virome_simulator.py simulator \
    -c step2_simulator/eval2a_config.txt \
    -o step2_simulator/eval2a_abundance --mode PE -l 150 \
    --seed 42 --threads 30

# ====== Step 4: 评估2B 已知检测(突变容忍) -- 5突变率 x 深度100x ======
for mut in 0 5 10 15 30; do
    mut_dir="step2_simulator/eval2b_mut_${mut}"

    # 1. 突变
    python sim/virome_simulator.py mutate \
        -i step1_eval_viruses/ -o "${mut_dir}" -r "${mut}"

    # 2. 配置
    python sim/virome_simulator.py gen-config \
        -i "${mut_dir}" \
        -o "step2_simulator/eval2b_mut${mut}_config.txt" \
        -d 100 -l 150

    # 3. 模拟测序
    python sim/virome_simulator.py simulator \
        -c "step2_simulator/eval2b_mut${mut}_config.txt" \
        -o "step2_simulator/eval2b_mut${mut}" \
        --mode PE -l 150 --seed 42 --threads 30
done

# 金标准生成
python deps/prepare_cami_input.py \
    -g step2_simulator/eval_known_virus_gold.tsv \
    -r final.cluster.ref_info.tsv \
    -o step2_simulator/eval_merged_input.tsv --pct
```

### 2.2.3.2 病毒鉴定
```bash
mkdir -p step6_sim && cd step6_sim
ln -s ../step2_simulator/eval2*gz ./
cd ..

# plant-virus数据库
PVDB=~/database/virus-db/plantvirus-db/plant_virus_db/5.virus.ref.build.db

bash eval/run_eval_known_virus.sh \
    --ref-fasta final.cluster.ref.fasta \
    --ref-info final.cluster.ref_info.tsv \
    --genes-cov virus_genes_cov.tsv \
    --ref-db $PVDB --multitool-db $PVDB \
    --threads 1 --align-threads 30 --multitool-threads 30 --multitool-jobs 6 \
    --input-dir step6_sim/ \
    --output-dir step6_identify_PVDB \
    --logdir step6_logs_PVDB

# 结果标准化
cd step6_identify_PVDB
for tool in kallisto salmon bowtie2 bwa minimap2 strobealign bwa-mem2 hisat2; do
    echo "Processing: $tool ..."
    python ../metrics/map2cami.py \
        -i $tool/summary/all_viruses.summary.tsv \
        --tool "$tool" \
        -o "${tool}.cami_profiles" \
        -v Asm_EM_Reads
done

# 合并 kmer 分类结果
for i in */; do
    d="${i%/}"
    python ../deps/batch_class.merge.py \
        --kraken2 $d/$d.kraken2.out \
        --kaiju $d/$d.kaiju.out \
        --kraken2x $d/$d.kraken2-x.out \
        --krakenuniq $d/$d.krakenuniq.out \
        --kunpeng $d/$d.kunpeng.output_1.txt \
        --centrifuger $d/$d.centrifuger.out \
        --metabuli $d/$d.metabuli_classifications.tsv \
        --sylph $d/$d.sylph.tsv \
        --sylph_map $PVDB/seqid2taxid.map \
        -tax $PVDB/ref.virus.build_taxonomy.txt \
        -outdir $d.merge -prefix $d --mode all
done

# 批量转换 混合报告
python ../deps/report2cami.py \
    -i multitool/*merge/*_tree_report.tsv \
    -o kmer.cami_profiles
cd ..
```

### 2.2.3.3 结果分析
```bash
mkdir -p step6_result_analysis && cd step6_result_analysis
mkdir -p ALL.cami_profiles
cp ../step6_identify_PVDB/*cami_profiles/* ALL.cami_profiles/

# 金标准生成
python ../metrics/make_cami.gold.py \
    -i ../step2_simulator/eval_merged_input.tsv \
    -s 0 -t 2 -r 3 -p 4 \
    -o ../step6_identify_PVDB/eval_cami_gold_results

# OPAL 评估
python ../metrics/run_opal_auto.py \
    --input-dir ALL.cami_profiles \
    --gold-dir ../step6_identify_PVDB/eval_cami_gold_results \
    --out-dir OPAL_Results_ALL_Datasets

# 箱线图
python ../plots/plot.cami.box.py \
    --input OPAL_Results_ALL_Datasets/results.tsv \
    --output Performance_Boxplots.strain.png --rank strain

# 雷达图
python ../plots/plot_opal_radar.py \
    -i OPAL_Results_ALL_Datasets/by_rank/ -f all \
    --out-dir OPAL_Results_ALL_Datasets/radar_plots

# 资源消耗图
python ../plots/plot_identify.resources.py \
    -i ../step6_identify_PVDB -o ../step6_identify_analysis/
cd ..
```

## 2.2.4 病毒组装基准分析

> 参考: Gupta & Kumar (2022), Forouzan et al. (2018), Meyer et al. (2022)
> 单病毒组装工具: haploflow, IVA, SAVAGE

### 2.2.4.2 数据模拟
```bash
# ====== Stage 1: 深度梯度 (无突变, 无背景) ======
VIRUS_DIR="step1_eval_viruses"

for depth in 1 5 10 20 50 100; do
    for rep in 1 2 3; do
        seed=$((42 + rep - 1))
        tag="depth${depth}_rep${rep}"
        [ -f "step7_sim/eval3_${tag}_PE_R1.fastq.gz" ] && continue

        python sim/virome_simulator.py gen-config \
            -i "$VIRUS_DIR" -d $depth -l 150 \
            -o "step2_simulator/eval3_config_${tag}.txt"

        python sim/virome_simulator.py simulator \
            -c "step2_simulator/eval3_config_${tag}.txt" \
            -o "step7_sim/eval3_${tag}" --mode PE -l 150 \
            --seed $seed --threads 30
    done
done
# 验证: 6 depths x 3 reps = 18
ls step7_sim/eval3_depth*_PE_R1.fastq.gz | wc -l

# ====== Stage 2: 宿主背景干扰 (50x深度, 0%突变, 6个背景比) ======
FIXED_DEPTH=50
VIRUS_DIR="step1_eval_viruses"

VIRUS_TOTAL_READS=$(python -c "
import os, math
total = 0
for f in os.listdir('$VIRUS_DIR'):
    if f.endswith('.fasta'):
        with open(os.path.join('$VIRUS_DIR', f)) as fh:
            for line in fh:
                if line.startswith('>'): continue
                total += len(line.strip())
reads = math.ceil(total * $FIXED_DEPTH / 150)
print(reads)
")
echo "Virus reads @ ${FIXED_DEPTH}x: $VIRUS_TOTAL_READS"

for ratio in 0.1 0.5 1 2 5 10; do
    bg_reads=$(python -c "print(int($VIRUS_TOTAL_READS * $ratio))")
    for rep in 1 2 3; do
        seed=$((42 + rep - 1))
        tag="bg_r${ratio}_rep${rep}"
        [ -f "step7_sim/eval3_${tag}_PE_R1.fastq.gz" ] && continue

        python sim/virome_simulator.py LoD_mix \
            -i "$VIRUS_DIR" --bgref ningxia.genome.fasta \
            --bg-reads $bg_reads --depths $FIXED_DEPTH \
            --mode PE -l 150 --seed $seed --threads 30 \
            -o "step7_sim/eval3_${tag}_tmp" --all-in-one

        mv "step7_sim/eval3_${tag}_tmp/Host_Depletion_Mixed_PE_R1.fastq.gz" \
           "step7_sim/eval3_${tag}_PE_R1.fastq.gz" 2>/dev/null
        mv "step7_sim/eval3_${tag}_tmp/Host_Depletion_Mixed_PE_R2.fastq.gz" \
           "step7_sim/eval3_${tag}_PE_R2.fastq.gz" 2>/dev/null
        rm -rf "step7_sim/eval3_${tag}_tmp"
    done
done
# 验证: 6 ratios x 3 reps = 18
ls step7_sim/eval3_bg*_PE_R1.fastq.gz | wc -l
echo "Total samples: $(ls step7_sim/eval3_*_PE_R1.fastq.gz 2>/dev/null | wc -l)"

# 合并参考FASTA
mkdir -p step7_ref
cat step1_eval_viruses/*.fasta > step7_ref/eval3_reference.fasta
echo "Reference: $(grep -c '^>' step7_ref/eval3_reference.fasta) genomes"
```

### 2.2.4.3 评估运行
```bash
bash eval/run_eval_assembly.sh \
    --sim-data step7_sim/ --output-dir step7_assembly \
    --batch 4 --threads 30 --jobs 4 --merge-jobs 5
```

### 2.2.4.4 基准测试与绘图
```bash
# 按模式分目录
mkdir -p step7_assembly_depth step7_assembly_bg
mv step7_assembly/eval3_depth* step7_assembly_depth/
mv step7_assembly/eval3_bg*   step7_assembly_bg/

for mode in bg depth; do
    python metrics/benchmark_mq.py \
        -i step7_assembly_${mode} -r step7_ref \
        -o step7_${mode}_benchmark --mode 6 -j 10 -t 20

    python metrics/benchmark_chimeric.py \
        -i step7_assembly_${mode} -r step7_ref \
        -o step7_${mode}_benchmark --mode 6 -j 10 -t 20

    python metrics/benchmark_resource.py \
        -i step7_assembly_${mode} \
        -o step7_${mode}_benchmark --mode 6

    python metrics/benchmark_summarize.py \
        -d step7_${mode}_benchmark --mode 6 --phase all \
        --viruses NC_002030.1 OR489165.1 AF395872
done
```

## 2.2.5 病毒鉴定基准分析

| 工具 | 原理 | 版本 |
|------|------|------|
| BLASTn/Diamond | 核酸/蛋白同源搜索 | E-value=1e-5, Top5 |
| VirBot | 蛋白HMM搜索 | v1.0 |
| VirSorter2 | 蛋白同源搜索 | v2.2.4 |
| ViralVerify | 蛋白HMM搜索 | v1.1 |
| RdRpCatch | RdRp结构域搜索 | v1.0 |
| ViraLM | Transformer语言模型 | 2025.01 |
| VirHunter | 深度学习 | generalistic |
| geNomad | 深度学习 | v1.9.0 |
| MetaBuli | K-mer频率 | v1.0 + RVDB-v31 |
| Ensemble | 多工具集成 | 投票/共识 |

### 2.2.5.1 数据模拟
```bash
# Pfam 注释 + 陷阱序列
python prep/prep_extract_conserved_traps.py \
    --host-proteins ningxia.pep.fasta \
    --host-genome ningxia.genome.fasta \
    --host-gff ningxia.genome.gff3 \
    --pfam-db ~/database/pfam-v35/Pfam-A.hmm \
    --outdir step3_conserved_traps/ --n-sequences 1000 --seed 42

# 鉴定评估测试集 (混合新病毒模拟)
python prep/prep_master_eval_dataset.py \
    --known-fasta step1_eval_viruses/ \
    --novel-fasta verified_new_plant_viruses.fasta \
    --extra-novel-fasta new.virus.fasta \
    --conserved-fasta step3_conserved_traps/ \
    --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
    --host-fasta ningxia.genome.fasta \
    --virus-meta step1_eval_viruses/selected_viruses.tsv \
    --known-mutations 0 5 10 15 \
    --coverage-levels 100 90 70 50 30 10 \
    --similarity 70 --n-per-coverage 5 --n-neg 1000 \
    --outdir step3_master_eval/ --seed 42
```

### 2.2.5.2 病毒鉴定分析
```bash
python deps/virus_identification.py \
    -i step3_master_eval/evaluation_sequences.fasta \
    -o step8_result \
    --identify_tools all --blast_mode both \
    --blast_top_n 5 --blast_evalue 1e-5 \
    --db_dir ~/database/virus-db \
    --virus_protein_db ~/database/virus-db/RVDB-v31/RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd \
    --uniprot_db ~/database/uniport_db/uniref90/uniref90.dmnd \
    --viroids_db ~/database/virus-db/viroids-db/viroids.fasta.blast.db \
    --virsorter_db ~/database/virus-db/virsorter2_db/ \
    --viralverify_hmm ~/database/virus-db/viralverify_db/nbc_hmms.hmm \
    --virhunter_path ~/biosoft/virus/virhunter/virhunter/predict_cpu.py \
    --virhunter_weights ~/biosoft/virus/virhunter/weights/generalistic \
    --metabuli_db ~/database/virus-db/RVDB-v31/RVDB_viroids.metabuli_db \
    --virbot_path ~/biosoft/virus/VirBot/VirBot.py \
    --viralm_path ~/bin/viralm_cpu.py \
    --virus_taxid ~/database/virus-db/taxIDs/viral_taxIDs.txt \
    --virsorter_group "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" \
    -j 10 -t 30

# 绘图分析
python metrics/eval_identification.py \
    --result-dir step8_result/step3_master_eval/ \
    --outdir step8_result_analysis/ \
    --labels step3_master_eval/sequence_labels_category.tsv \
    --virus-dir step1_eval_viruses/ \
    --min-virus-length 500 --filter-mode all \
    --prefix step3_master_eval_virus
```

## 2.2.6 病毒分类基准分析

### 2.2.6.1 数据模拟
```bash
# 分类评估测试集 (一次性生成完整数据集)
python deps/prep_build_class_eval_seqs.py \
    --virus-dir step1_eval_viruses/ \
    --ref-info final.cluster.ref_info.tsv \
    --ref-fasta final.cluster.ref.fasta \
    --novel-virus novel_viruses.2026.fasta \
    --conserved-fasta step3_conserved_traps/ \
    --eve-fasta C-RVDBvCurrent.PLN.not-virus.fasta \
    --host-fasta ningxia.genome.fasta \
    --n-decoys 300 \
    --coverage-levels 100 90 80 70 60 50 40 --n-per-cov 5 \
    --mutation-rates 0.00 0.05 0.10 0.15 --n-per-mut 1 \
    --email 1771182368@qq.com \
    --outdir step4_classification_eval/ --seed 42

# 补充 NCBI 元数据
python metrics/enrich_metadata_taxonomy.py \
    --meta step4_classification_eval/test_metadata.tsv \
    --out step4_classification_eval/test_metadata_full.tsv \
    --email 1771182368@qq.com
```

> 预期输出: 18586 test sequences
> - known: 7104 seqs (mut0/5/10/15, cov40-100)
> - novel: 10582 seqs (mut0/5/10/15, cov40-100)
> - pfam/eve/host: 各300 seqs (decoy)

### 2.2.6.2 分类测试
```bash
python deps/virus_classifier.py \
    --genomes step4_classification_eval/evaluation_sequences.fasta \
    -s class \
    --output-dir step9_classification/classifier \
    --tools all --threads 20 \
    --uniprot-db ~/database/virus-db/RVDB-30/U-RVDBv30_prot.dmnd \
    --genomad-db ~/database/virus-db/genomad_db/ \
    --metabuli-db ~/database/virus-db/RVDB-30/5.virus.ref.build.db/ref.virus.build.metabuli_db/ \
    --cat-db ~/database/virus-db/RVDB-30/CAT-db/db/ \
    --cat-tax ~/database/virus-db/RVDB-30/CAT-db/tax/ \
    --mmseqs-db ~/database/virus-db/RVDB-30/RVDB.mmseqs \
    --vitap-db ~/database/virus-db/vitap-db/VMR-MSL40_DB \
    --acvirus-db ~/database/virus-db/acvirus_db \
    --vcontact3-db ~/database/virus-db/vConTACT3_db

# 合并分析
Rscript deps/virus_classifier_analysis14.R \
    --combined step9_classification/classifier/classifier_combined_taxonomy.tsv \
    -o step9_classification/integrated/
```

### 2.2.6.3 评估
```bash
# 合并分类评估
python metrics/run_full_analysis.py \
    --predictions-dir step9_classification/integrated/ \
    --meta step4_classification_eval/test_metadata_full.tsv \
    --outdir step9_classification/analysis/

# 工具组合优化
python metrics/opt_classifier_ensemble.py \
    --predictions step9_classification/integrated/ \
    --meta step4_classification_eval/test_metadata_full.tsv \
    --outdir step9_classification/ensemble_opt

python metrics/ensemble_summarize.py \
    --opt-tsv step9_classification/ensemble_opt/ensemble_optimization.tsv \
    --outdir step9_classification/ensemble_opt/

# 绘图
python plots/plot_classification_comparison.py \
    --input step9_classification/analysis/ \
    --outdir step9_classification/analysis/plots/

python metrics/eval_segmented_virus.py \
    --selected step1_eval_viruses/selected_viruses.tsv \
    --meta step4_classification_eval/test_metadata_full.tsv \
    --outdir step9_classification/analysis/segmented/
```
