#!/usr/bin/env python3
"""
eval_plant_virus 集中配置文件
===============================
所有路径、数据库、外部工具的唯一定义点。
可通过环境变量 EVAL_<KEY> 覆盖任何路径。

用法:
    from config import DEPS_DIR, DATABASE_PATHS, ...
    python config.py --export-sh   # 导出为 shell 可 source 的格式
"""

import os

# ── 项目根 ────────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# ── 外部依赖目录 (eval_plant_virus/deps/bin/) ─────────────
DEPS_DIR = os.path.join(PROJECT_ROOT, "deps")

def _d(name):
    return os.path.join(DEPS_DIR, name)

# 核心管线脚本
ASSEMBLY_PIPELINE           = _d("assembly_pipeline.py")
HOST_DEPLETION              = _d("host_depletion.py")
VIRUS_CLASSIFIER            = _d("virus_classifier.py")
VIRUS_IDENTIFICATION        = _d("virus_identification.py")
BATCH_CLASS_MERGE           = _d("batch_class.merge.py")
BATCH_CLASS_READS           = _d("batch_class.reads.py")
BATCH_VIRUS_DEPTH           = _d("batch_virus_depth40.py")
PREPARE_CAMI_INPUT          = _d("prepare_cami_input.py")
REPORT2CAMI                 = _d("report2cami.py")
KRAKEN2CAMI                 = _d("kraken2cami.py")
RUN_CLUSTER_PIPELINE        = _d("run_cluster_pipeline.py")
PREP_BUILD_CLASS_EVAL_SEQS  = _d("prep_build_class_eval_seqs.py")
VIRUS_CLASSIFIER_ANALYSIS_R = _d("virus_classifier_analysis14.R")
PLOT_CAMI_BOX               = os.path.join(PROJECT_ROOT, "plots", "plot.cami.box.py")

# 项目内脚本
DETECT_CHIMERIC_CONTIGS = os.path.join(PROJECT_ROOT, "detect_chimeric_contigs.py")
VIROME_SIMULATOR        = os.path.join(PROJECT_ROOT, "sim", "virome_simulator.py")

# ── 数据库根路径 ────────────────────────────────────────
DB_ROOT = os.path.expanduser("~/database")

DATABASE_PATHS = {
    # Plant virus DB
    "ref_info":     os.path.join(DB_ROOT, "virus-db", "plantvirus-db", "plant_virus_db",
                                 "3.final-ref-virus.db", "final.cluster.ref_info.tsv"),
    "ref_fasta":    os.path.join(DB_ROOT, "virus-db", "plantvirus-db", "plant_virus_db",
                                 "3.final-ref-virus.db", "final.cluster.ref.fasta"),
    "plant_virus_db": os.path.join(DB_ROOT, "virus-db", "plantvirus-db", "plant_virus_db",
                                   "5.virus.ref.build.db"),
    "genes_cov":    os.path.join(DB_ROOT, "virus-db", "plantvirus-db", "plant_virus_db",
                                 "5.virus.ref.build.db", "virus_genes_cov.tsv"),

    # RVDB
    "rvdb_prot":    os.path.join(DB_ROOT, "virus-db", "RVDB-v31",
                                 "RVDB_viroids.diamond_db", "U-RVDBv31.0-prot_unique.dmnd"),
    "rvdb_metabuli": os.path.join(DB_ROOT, "virus-db", "RVDB-v31", "RVDB_viroids.metabuli_db"),
    "rvdb_eve":     os.path.join(DB_ROOT, "virus-db", "RVDB-v31",
                                 "C-RVDBvCurrent.extract", "C-RVDBvCurrent.PLN.not-virus.fasta"),

    # Classification DBs
    "mmseqs_db":    os.path.join(DB_ROOT, "virus-db", "RVDB-30", "RVDB.mmseqs"),
    "vitap_db":     os.path.join(DB_ROOT, "virus-db", "vitap-db", "VMR-MSL40_DB"),
    "acvirus_db":   os.path.join(DB_ROOT, "virus-db", "acvirus_db"),
    "vcontact3_db": os.path.join(DB_ROOT, "virus-db", "vConTACT3_db"),
    "genomad_db":   os.path.join(DB_ROOT, "virus-db", "genomad_db"),
    "cat_db":       os.path.join(DB_ROOT, "virus-db", "RVDB-30", "CAT-db", "db"),
    "cat_tax":      os.path.join(DB_ROOT, "virus-db", "RVDB-30", "CAT-db", "tax"),

    # Other DBs
    "uniref90":     os.path.join(DB_ROOT, "uniport_db", "uniref90", "uniref90.dmnd"),
    "viroids_db":   os.path.join(DB_ROOT, "virus-db", "viroids-db", "viroids.fasta.blast.db"),
    "pfam":         os.path.join(DB_ROOT, "pfam-v35", "Pfam-A.hmm"),
    "virsorter2":   os.path.join(DB_ROOT, "virus-db", "virsorter2_db"),
    "viralverify":  os.path.join(DB_ROOT, "virus-db", "viralverify_db", "nbc_hmms.hmm"),
    "taxid":        os.path.join(DB_ROOT, "virus-db", "taxIDs", "viral_taxIDs.txt"),

    # Host
    "host_k2":      os.path.join(DB_ROOT, "host_db", "kraken2"),
    "host_hisat2":  os.path.join(DB_ROOT, "host_db", "hisat2", "host"),
    "host_genome":  os.path.join(DB_ROOT, "genome-data", "ningxia-genome", "1.genome",
                                 "ningxia.genome.fasta"),
    "host_pep":     os.path.join(DB_ROOT, "genome-data", "ningxia-genome", "1.genome",
                                 "ningxia.pep.fasta"),
    "host_gff":     os.path.join(DB_ROOT, "genome-data", "ningxia-genome", "1.genome",
                                 "ningxia.genome.gff3"),
}

# ── 生物软件路径 ──────────────────────────────────────────
BIOSOFT_DIR = os.path.expanduser("~/biosoft")

TOOL_PATHS = {
    "virbot":     os.path.join(BIOSOFT_DIR, "virus", "VirBot", "VirBot.py"),
    "virhunter":  os.path.join(BIOSOFT_DIR, "virus", "virhunter", "virhunter", "predict_cpu.py"),
    "virhunter_w": os.path.join(BIOSOFT_DIR, "virus", "virhunter", "weights", "generalistic"),
    "mutation_sim": os.path.join(BIOSOFT_DIR, "binary", "mutation-simulator"),
}

# Conda / 系统工具
MAMBAFORGE_BIN = os.path.expanduser("~/mambaforge/bin")
TOOL_SEARCH_PATHS = [
    os.path.join(BIOSOFT_DIR, "binary"),
    MAMBAFORGE_BIN,
    os.path.expanduser("~/bin"),
]


# ── 环境变量覆盖 ──────────────────────────────────────────
def get(key, default):
    """允许通过环境变量 EVAL_<KEY> 覆盖任何路径"""
    return os.environ.get(f"EVAL_{key.upper()}", default)


# ── Shell 导出 ────────────────────────────────────────────
def export_sh():
    """导出可在 shell 中 source 的配置"""
    lines = [
        f'export EVAL_ROOT="{PROJECT_ROOT}"',
        f'export EVAL_DEPS="{DEPS_DIR}"',
        f'export EVAL_DB_ROOT="{DB_ROOT}"',
    ]
    for k, v in DATABASE_PATHS.items():
        lines.append(f'export EVAL_DB_{k.upper()}="{v}"')
    for k, v in TOOL_PATHS.items():
        lines.append(f'export EVAL_TOOL_{k.upper()}="{v}"')
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--export-sh", action="store_true", help="导出 shell 格式")
    p.add_argument("--check", action="store_true", help="检查所有路径存在性")
    args = p.parse_args()

    if args.export_sh:
        print(export_sh())
    elif args.check:
        missing = []
        for k, v in DATABASE_PATHS.items():
            if not os.path.exists(v):
                missing.append(f"  DB  {k}: {v}")
        for k, v in TOOL_PATHS.items():
            if not os.path.exists(v):
                missing.append(f"  TOOL {k}: {v}")
        for name in os.listdir(DEPS_DIR) if os.path.isdir(DEPS_DIR) else []:
            pass  # deps/bin may not exist yet
        if missing:
            print("缺少的路径:")
            for m in missing:
                print(m)
        else:
            print("所有路径检查通过。")
