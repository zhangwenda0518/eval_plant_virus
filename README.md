# eval_plant_virus — 植物病毒检测流程基准评估工具集

植物病毒生物信息学检测流程的标准化基准评估体系（第一章配套脚本）。

## 目录结构

```
eval_plant_virus/
├── README.md                           # 本文件
├── config/default_paths.sh             # 默认数据库路径配置
├── prep/                               # 数据准备脚本
│   ├── prep_select_eval_viruses.py     #  选取评估用病毒（非节段RefSeq完整基因组）
│   ├── prep_build_id_eval_seqs.py      #  构建鉴定评估数据集（覆盖度梯度+EVE/陷阱）
│   ├── prep_build_class_eval_seqs.py   #  构建分类评估数据集（防信息泄漏）
│   └── prep_extract_conserved_traps.py #  提取保守结构域陷阱序列
├── sim/                                # 模拟测序脚本
│   └── virome_simulator.py             #  全功能模拟数据生成引擎
├── eval/                               # 评估调度脚本
│   ├── run_eval_known_virus.sh         #  评估一：已知病毒检测方法比较
│   ├── run_eval_assembly.sh            #  评估二：病毒组装方法比较
│   ├── run_eval_host_depletion.sh      #  宿主过滤消融实验
│   ├── run_eval_identification.sh      #  评估三：候选病毒鉴定策略比较
│   └── run_eval_classification.sh      #  评估四：病毒分类方法比较
├── metrics/                            # 评估指标计算脚本
│   ├── detect_chimeric_contigs.py      #  嵌合Contig检测
│   ├── calc_auprc.py                   #  PR曲线/AUPRC计算
│   ├── calc_abundance_accuracy.py      #  丰度定量准确性(Bray-Curtis/RMSE/Spearman)
│   └── calc_classification_stratified.py # 分类按覆盖度分层报告
└── report/                             # 汇总报告脚本（待完成）
    └── summarize_benchmark.py
```

## 快速开始

### 前置依赖

```bash
conda install -c bioconda art iss seqkit bbtools pandas numpy biopython scipy scikit-learn
pip install matplotlib seaborn
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

## 四个评估环节

| 评估 | 问题 | 主要指标 | 统计方法 |
|------|------|---------|---------|
| 评估一：已知病毒检测 | 比对 vs 准映射 vs k-mer 孰优孰劣？ | F1/LoD/丰度准确性(Bray-Curtis,RMSE,Spearman) | Friedman+Nemenyi |
| 评估二：病毒组装 | MEGAHIT vs SPAdes vs Penguin，分开 vs 合并？ | Genome fraction/NGA50/嵌合率/完全组装比 | Wilcoxon |
| 评估三：候选病毒鉴定 | 三种搜索原理互补性？对抗策略效果？ | AUPRC/Precision/Recall/F1 | — |
| 评估四：病毒分类 | MMseqs2 vs VITAP vs ACVirus？整合增益？ | 分覆盖度分层科/属/种准确率 | Cochran's Q+McNemar |

## 附加评估

- **宿主过滤消融实验**：Kraken2/HISAT2/RiboDetector 各步骤独立贡献
- **嵌合体率检测**：组装质量交叉验证指标
- **覆盖度梯度分层**：分类工具在不同序列完整度下的表现
