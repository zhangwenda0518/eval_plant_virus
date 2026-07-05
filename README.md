# eval_plant_virus — 植物病毒检测流程基准评估工具集

植物病毒生物信息学检测流程的标准化基准评估体系（第一章配套脚本）。

## 目录结构

```
eval_plant_virus/
├── README.md                           # 本文件
├── config.py                           # 集中路径配置 (替代旧的 config/default_paths.sh)
├── environment.yml                     # Conda 环境依赖
├── deps/bin/                           # 外部管线依赖 (13个脚本, 从 ~/bin/ 归入)
│   ├── assembly_pipeline.py
│   ├── host_depletion.py
│   ├── virus_classifier.py             # <- MMPV-RNA 最新版
│   ├── virus_identification.py         # <- MMPV-RNA 最新版
│   ├── batch_class.merge.py / batch_class.reads.py
│   ├── batch_virus_depth40.py
│   ├── prepare_cami_input.py / report2cami.py / kraken2cami.py
│   ├── run_cluster_pipeline.py
│   └── prep_build_class_eval_seqs.py / virus_classifier_analysis14.R
├── prep/                               # 数据准备脚本
│   ├── prep_select_eval_viruses.py     #  选取评估用病毒
│   ├── prep_build_id_eval_seqs.py      #  构建鉴定评估数据集
│   ├── prep_build_class_eval_seqs.py   #  构建分类评估数据集
│   ├── prep_extract_conserved_traps.py #  保守结构域陷阱序列
│   ├── prep_master_eval_dataset.py     #  主评估数据集生成
│   ├── prep_master_eval.py             #  主评估入口
│   └── prep_simulation_datasets.py     #  模拟数据集
├── sim/                                # 模拟测序
│   └── virome_simulator.py
├── eval/                               # 评估调度
│   ├── run_eval_known_virus.sh         #  评估一：已知病毒检测
│   ├── run_eval_assembly.sh            #  评估二：病毒组装
│   ├── run_eval_host_depletion.sh      #  宿主过滤消融
│   ├── run_eval_identification.sh      #  评估三：候选病毒鉴定
│   └── run_eval_classification.sh      #  评估四：病毒分类
├── metrics/                            # 评估指标
│   ├── benchmark_mq.py / benchmark_chimeric.py / benchmark_resource.py / benchmark_summarize.py
│   ├── eval_identification.py / eval_identification2.py
│   ├── eval_host_depletion.py / eval_host_prediction.py
│   ├── eval_cluster_pipeline.py / eval_segmented_virus.py
│   ├── run_full_analysis.py / opt_classifier_ensemble.py
│   └── calc_auprc.py / calc_abundance_accuracy.py / calc_classification_stratified.py
└── plots/                              # 绘图
    ├── plot.cami.box.py / plot_opal_radar.py
    ├── plot_host_depletion.py / plot_host_depletion.resources.py
    ├── plot_classification_comparison.py
    └── plot_identification_comparison.py / plot_identify.resources.py
```

## 快速开始

### 前置依赖

```bash
# 一键创建 conda 环境
conda env create -f environment.yml
conda activate eval_plant_virus

# 检查路径配置
python config.py --check
```

### Step 1: 选取50个评估用病毒

```bash
python prep/prep_select_eval_viruses.py \
    --ref-info ~/database/virus-db/plant_virus_db/final.cluster.ref_info.tsv \
    --ref-fasta ~/database/virus-db/plant_virus_db/final.cluster.ref.fasta \
    --n-viruses 50 --outdir eval_viruses_50/ --seed 42
```

### Step 2: 生成模拟数据

```bash
python sim/virome_simulator.py benchmark \
    -i eval_viruses_50/ --bgref host_genome.fasta \
    --depth 1000 --mut-rates 0 15 --depths 200k 500k 1M 2M --repeats 5 \
    --mode PE -l 150 --seed 42 --threads 40 -o sim_data/
```

### Step 3: 构建评估测试集

```bash
# 鉴定评估集（含EVE/转座子负样本）
python prep/prep_build_id_eval_seqs.py \
    --virus-fasta eval_viruses_50/ --coverage-levels 100 80 60 40 20 \
    --n-per-coverage 2 --host-genome host_genome.fasta \
    --conserved-prots conserved_traps/ \
    --eve-fasta host_eves.fasta \
    --outdir eval_identification/

# 分类评估集
python prep/prep_build_class_eval_seqs.py \
    --virus-dir eval_viruses_50/ --ref-info ref_info.tsv --ref-fasta ref.fasta \
    --coverage-levels 100 80 60 40 20 --n-per-coverage 2 \
    --outdir eval_classification/
```

### Step 4-7: 运行四个评估

```bash
bash eval/run_eval_known_virus.sh
bash eval/run_eval_host_depletion.sh
bash eval/run_eval_assembly.sh
bash eval/run_eval_identification.sh
bash eval/run_eval_classification.sh
```

### Step 8: 计算补充指标

```bash
# 嵌合体检测
python metrics/detect_chimeric_contigs.py -i contigs.blastn.tsv -o chimeras.tsv

# AUPRC
python metrics/calc_auprc.py --predictions results.tsv --labels labels.tsv --out auprc/

# 丰度准确性
python metrics/calc_abundance_accuracy.py --predictions quant.tsv --gold gold.tsv --out accuracy/

# 分类分层报告
python metrics/calc_classification_stratified.py --predictions integrated.tsv --meta test_meta.tsv --out stratified/
```

## 七个评估环节

| 评估 | 调度脚本 | 问题 | 主要指标 |
|------|---------|------|---------|
| 评估一：宿主过滤消融 | `run_eval_host_depletion.sh` | Kraken2/HISAT2 各步骤独立贡献？ | 病毒保留率/宿主去除率/资源消耗 |
| 评估二：已知病毒检测 | `run_eval_known_virus.sh` | 比对 vs 准映射 vs k-mer 孰优？ | F1/LoD/丰度准确性(Bray-Curtis,RMSE,Spearman) |
| 评估三：病毒组装 | `run_eval_assembly.sh` | MEGAHIT vs SPAdes vs Penguin？ | Genome fraction/NGA50/嵌合率 |
| 评估四：候选病毒鉴定 | `run_eval_identification.sh` | 多原理互补性？对抗策略效果？ | AUPRC/Precision/Recall/F1 |
| 评估五：病毒分类 | `run_eval_classification.sh` | MMseqs2 vs VITAP vs ACVirus？ | 分覆盖度分层科/属/种准确率 |
| 评估六：序列去重聚类 | `run_dedup_clustering.sh` | mmseqs vs vclust 多种模式？ | ARI/V-measure/NMI |
| 评估七：宿主分类基准 | `run_eval_host_prediction.sh` | RNAVirHost/PhaBOX2/ICTV 集成？ | Precision/Recall/F1/Accuracy |
