#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
  dereplication.py — 统一的病毒基因组去冗余评估模块

  整合两种去冗余方法:
    1. BLAST+aniclust  — 基于 BLASTn 全对比对 + 加权 ANI + 贪婪质心聚类
    2. dRep             — MASH预筛选 + ANIm精确聚类 (带病毒适应参数)

  用法:
    python dereplication.py -i input.fasta -o outdir -m both -t 8

  依赖:
    - BLAST+  (makeblastdb, blastn)
    - dRep    (可选, 仅 --method drep/both 时需要)
    - Python >= 3.7, BioPython
=============================================================================
"""

import os, sys, gzip, time, json, logging, argparse, shutil, subprocess, tempfile
from dataclasses import dataclass, field, asdict
from collections import defaultdict
from typing import Dict, List, Set, Optional, Tuple
from pathlib import Path

# ============================================================================
# 病毒适应的 dRep 参数 (从 ViOTUcluster & ViWrap 提取)
# ============================================================================
DREP_VIRUS_PARAMS = {
    "pa": 0.8,          # primary ANI (MASH预筛选)
    "sa": 0.95,         # secondary ANI (种水平阈值 ≈ ICTV)
    "nc": 0.85,         # min overlap 85% = MIUViG标准
    "ignore_genome_quality": True,  # 病毒基因通常不完整, CheckM无意义
    "comW": 0,          # 完整度权重 = 0
    "conW": 0,          # 污染度权重 = 0
    "strW": 0,          # 菌株异质性 = 0 (病毒不适用)
    "N50W": 0,          # N50 权重 = 0 (多片段病毒无意义)
    "sizeW": 1,         # 唯一权重: 基因组大小
    "centW": 0,         # 中心度 = 0
    "min_length": 3000, # 最小基因组长度
    "skip_plots": True, # 减少依赖
}

# ============================================================================
# 数据结构
# ============================================================================

@dataclass
class ClusteringResult:
    """聚类结果"""
    method: str                                    # "blast_aniclust" | "drep"
    clusters: Dict[str, List[str]] = field(default_factory=dict)  # rep → members
    representatives: List[str] = field(default_factory=list)
    n_clusters: int = 0
    n_singletons: int = 0
    n_total_sequences: int = 0
    mean_cluster_size: float = 0.0
    cluster_sizes: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("clusters", None)  # 不序列化完整cluster
        return d


@dataclass
class ComparisonReport:
    """两种方法的对比报告"""
    blast_result: Optional[ClusteringResult] = None
    drep_result: Optional[ClusteringResult] = None
    shared_representatives: Set[str] = field(default_factory=set)
    blast_only_representatives: Set[str] = field(default_factory=set)
    drep_only_representatives: Set[str] = field(default_factory=set)
    summary: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "blast": self.blast_result.to_dict() if self.blast_result else None,
            "drep": self.drep_result.to_dict() if self.drep_result else None,
            "n_shared": len(self.shared_representatives),
            "n_blast_only": len(self.blast_only_representatives),
            "n_drep_only": len(self.drep_only_representatives),
            "summary": self.summary,
        }


# ============================================================================
# 内嵌 anicalc.py 算法 (移植自 ViOTUcluster/anicalc.py)
# ============================================================================

def _parse_blast_line(line: str) -> dict:
    """解析 BLAST outfmt '6 std qlen slen' 的一行"""
    r = line.rstrip('\n').split('\t')
    return {
        "qname": r[0], "tname": r[1], "pid": float(r[2]),
        "len": float(r[3]),
        "qcoords": sorted([int(r[6]), int(r[7])]),
        "tcoords": sorted([int(r[8]), int(r[9])]),
        "qlen": float(r[-2]), "tlen": float(r[-1]),
        "evalue": float(r[-4]),
    }


def _yield_alignment_blocks(handle):
    """按 (qname, tname) 分组相邻的 HSP 行"""
    key, alns = None, None
    for line in handle:
        if line.startswith('#') or not line.strip():
            continue
        aln = _parse_blast_line(line)
        new_key = (aln["qname"], aln["tname"])
        if key is None:
            key, alns = new_key, [aln]
        elif new_key == key:
            alns.append(aln)
        else:
            if alns:
                yield alns
            key, alns = new_key, [aln]
    if alns:
        yield alns


def _prune_alns(alns: List[dict], min_length: int = 0, min_evalue: float = 1e-3) -> List[dict]:
    """剔除短HSP、高E-value、超出query全长110%的比对"""
    keep = []
    cur_aln = 0
    qry_len = alns[0]["qlen"]
    for aln in alns:
        qcoords = aln["qcoords"]
        aln_len = max(qcoords) - min(qcoords) + 1
        if aln_len < min_length or aln["evalue"] > min_evalue:
            continue
        if cur_aln >= qry_len or aln_len + cur_aln >= 1.10 * qry_len:
            break
        keep.append(aln)
        cur_aln += aln_len
    return keep


def _compute_ani(alns: List[dict]) -> float:
    """加权 ANI: Σ(len × pid) / Σ(len)"""
    return round(
        sum(a["len"] * a["pid"] for a in alns) / sum(a["len"] for a in alns), 2
    )


def _merge_coords(coords_list: List[List[int]]) -> int:
    """合并重叠坐标区间, 返回总比对长度"""
    coords = sorted(coords_list)
    nr = [coords[0]]
    for start, stop in coords[1:]:
        if start <= (nr[-1][1] + 1):
            nr[-1][1] = max(nr[-1][1], stop)
        else:
            nr.append([start, stop])
    return sum(stop - start + 1 for start, stop in nr)


def _compute_cov(alns: List[dict]) -> Tuple[float, float]:
    """计算 query覆盖率 和 target覆盖率 (去重合并HSP区间)"""
    qcov = round(100.0 * _merge_coords([a["qcoords"] for a in alns]) / alns[0]["qlen"], 2)
    tcov = round(100.0 * _merge_coords([a["tcoords"] for a in alns]) / alns[0]["tlen"], 2)
    return qcov, tcov


# ============================================================================
# 内嵌 aniclust.py 算法 (移植自 ViOTUcluster/aniclust.py)
# ============================================================================

def _parse_fasta_lengths(fasta_path: str, min_length: int = 1) -> Dict[str, int]:
    """读取FASTA, 返回 {seq_id: length}, 按长度降序"""
    seqs = {}
    try:
        from Bio import SeqIO
        for rec in SeqIO.parse(fasta_path, "fasta"):
            if len(rec.seq) >= min_length:
                seqs[rec.id] = len(rec.seq)
    except ImportError:
        # 纯Python fallback
        with open(fasta_path) as fh:
            sid, seq = None, []
            for line in fh:
                if line.startswith('>'):
                    if sid and len(''.join(seq)) >= min_length:
                        seqs[sid] = len(''.join(seq))
                    sid = line[1:].strip().split()[0]
                    seq = []
                else:
                    seq.append(line.strip())
            if sid and len(''.join(seq)) >= min_length:
                seqs[sid] = len(''.join(seq))
    return seqs


def _centroid_cluster(seq_lengths: Dict[str, int],
                      ani_file: str,
                      min_ani: float = 95,
                      min_tcov: float = 85,
                      min_qcov: float = 0) -> Dict[str, List[str]]:
    """
    贪婪质心聚类 (aniclust 算法)

    1. 按序列长度降序排列
    2. 最长未分配的序列成为新簇质心
    3. 所有满足阈值 (ANI/qcov/tcov) 的邻居归入当前簇
    4. 已分配序列不再参与后续质心竞争
    """
    # 1. 排序
    sorted_seqs = sorted(seq_lengths.keys(), key=lambda x: seq_lengths[x], reverse=True)
    logger = logging.getLogger(__name__)

    # 2. 读取ANI边
    edges: Dict[str, List[str]] = {s: [] for s in sorted_seqs}
    open_func = gzip.open if ani_file.endswith('.gz') else open
    with open_func(ani_file, 'rt') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            qname, tname, _, ani, qcov, tcov = parts[:6]
            try:
                ani_val, qcov_val, tcov_val = float(ani), float(qcov), float(tcov)
            except ValueError:
                continue
            if qname == tname:
                continue
            if qname not in edges or tname not in edges:
                continue
            if qcov_val < min_qcov or tcov_val < min_tcov or ani_val < min_ani:
                continue
            edges[qname].append(tname)

    # 3. 聚类
    clust_to_seqs: Dict[str, List[str]] = {}
    seq_to_clust: Dict[str, str] = {}

    for seq_id in sorted_seqs:
        if seq_id in seq_to_clust:
            continue
        # 新簇质心
        clust_to_seqs[seq_id] = [seq_id]
        seq_to_clust[seq_id] = seq_id
        for mem_id in edges[seq_id]:
            if mem_id not in seq_to_clust:
                clust_to_seqs[seq_id].append(mem_id)
                seq_to_clust[mem_id] = seq_id

    logger.info(f"  aniclust: {len(clust_to_seqs)} clusters from {len(sorted_seqs)} sequences")
    return clust_to_seqs


# ============================================================================
# FASTA 工具
# ============================================================================

def _write_representatives(fasta_path: str, representatives: Set[str], output_path: str):
    """从输入FASTA中提取代表序列写入新文件"""
    count = 0
    try:
        from Bio import SeqIO
        with open(output_path, 'w') as out:
            for rec in SeqIO.parse(fasta_path, "fasta"):
                if rec.id in representatives:
                    SeqIO.write(rec, out, "fasta")
                    count += 1
    except ImportError:
        write_flag = False
        with open(fasta_path) as fh, open(output_path, 'w') as out:
            for line in fh:
                if line.startswith('>'):
                    sid = line[1:].strip().split()[0]
                    write_flag = sid in representatives
                if write_flag:
                    out.write(line)
                    if line.startswith('>'):
                        count += 1
    return count


def _parse_fasta_headers(fasta_path: str) -> List[str]:
    """提取FASTA中所有序列ID"""
    ids = []
    try:
        from Bio import SeqIO
        for rec in SeqIO.parse(fasta_path, "fasta"):
            ids.append(rec.id)
    except ImportError:
        with open(fasta_path) as fh:
            for line in fh:
                if line.startswith('>'):
                    ids.append(line[1:].strip().split()[0])
    return ids


# ============================================================================
# dRep 结果解析
# ============================================================================

def _parse_drep_clusters(drep_outdir: str, min_length: int) -> ClusteringResult:
    """解析 dRep 输出 (Cdb.csv + Wdb.csv + Bdb.csv) 构建 ClusteringResult"""
    logger = logging.getLogger(__name__)
    result = ClusteringResult(method="drep")

    # 逐个子目录解析
    data_dirs = list(Path(drep_outdir).glob("*/data_tables"))
    if not data_dirs:
        logger.warning("  dRep: no data_tables found in output")
        return result

    for dt_dir in data_dirs:
        cdb_path = dt_dir / "Cdb.csv"
        wdb_path = dt_dir / "Wdb.csv"

        if not cdb_path.exists() or not wdb_path.exists():
            continue

        # Wdb: cluster → winning genome (代表)
        cluster_to_rep = {}
        with open(wdb_path) as fh:
            for line in fh:
                if line.startswith('genome'):
                    continue
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    rep = parts[0].rsplit('.', 1)[0]  # 去后缀
                    cluster_to_rep[parts[1]] = rep

        # Cdb: genome → cluster
        gn_to_cluster = {}
        with open(cdb_path) as fh:
            for line in fh:
                if line.startswith('genome'):
                    continue
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    gn = parts[0].rsplit('.', 1)[0]  # 去后缀
                    gn_to_cluster[gn] = parts[1]

        # 构建 cluster → members
        for cluster_id, rep in cluster_to_rep.items():
            members = [gn for gn, cl in gn_to_cluster.items() if cl == cluster_id]
            result.clusters[rep] = members
            result.cluster_sizes[rep] = len(members)

    result.representatives = list(result.clusters.keys())
    result.n_clusters = len(result.clusters)
    result.n_singletons = sum(1 for v in result.clusters.values() if len(v) == 1)
    result.n_total_sequences = sum(len(v) for v in result.clusters.values())
    result.mean_cluster_size = (result.n_total_sequences / result.n_clusters
                                if result.n_clusters > 0 else 0)

    logger.info(f"  dRep result: {result.n_clusters} clusters, "
                f"{result.n_singletons} singletons, "
                f"{result.n_total_sequences} total sequences")
    return result


# ============================================================================
# 主类: Dereplicator
# ============================================================================

class Dereplicator:
    """
    病毒基因组去冗余器

    Parameters
    ----------
    work_dir : str
        工作目录 (用于存放临时文件)
    threads : int
        线程数
    min_length : int
        最小序列/基因组长度
    min_ani : float
        ANI阈值 (默认95%)
    min_tcov : float
        最小target覆盖率 (默认85%)
    min_qcov : float
        最小query覆盖率 (默认0%)
    """

    def __init__(self, work_dir: str = "./dereplication_tmp",
                 threads: int = 8,
                 min_length: int = 3000,
                 min_ani: float = 95.0,
                 min_tcov: float = 85.0,
                 min_qcov: float = 0.0):
        self.work_dir = Path(os.path.abspath(work_dir))
        self.threads = threads
        self.min_length = min_length
        self.min_ani = min_ani
        self.min_tcov = min_tcov
        self.min_qcov = min_qcov
        self.logger = logging.getLogger(__name__)
        os.makedirs(self.work_dir, exist_ok=True)
        # 检查BLAST可用性
        self._blast_available = shutil.which("makeblastdb") is not None and shutil.which("blastn") is not None
        self._drep_available = shutil.which("dRep") is not None
        if not self._blast_available:
            self.logger.warning("BLAST+ (makeblastdb/blastn) not found; BLAST+aniclust method will fail.")
        if not self._drep_available:
            self.logger.warning("dRep not found; dRep method will fail.")

    # ----------------------------------------------------------------
    # 方法1: BLAST + anicalc + aniclust
    # ----------------------------------------------------------------

    def run_blast_aniclust(self, input_fasta: str, output_dir: str) -> ClusteringResult:
        """
        运行 BLAST+aniclust 去冗余流程

        Steps:
          makeblastdb → blastn (全vs全) → ANI计算 → 质心聚类 → 提取代表序列

        Parameters
        ----------
        input_fasta : str
            输入FASTA文件路径
        output_dir : str
            输出目录

        Returns
        -------
        ClusteringResult
        """
        os.makedirs(output_dir, exist_ok=True)
        fasta_name = Path(input_fasta).stem
        tmp_dir = os.path.join(self.work_dir, "blast_work")
        os.makedirs(tmp_dir, exist_ok=True)

        # BLAST工具在非ASCII路径下可能失败 (Windows), 回退到系统临时目录
        try:
            tmp_dir.encode('ascii')
        except UnicodeEncodeError:
            tmp_dir = tempfile.mkdtemp(prefix="derep_blast_")
            os.makedirs(tmp_dir, exist_ok=True)
            self.logger.info(f"  Using temp dir (ASCII-safe): {tmp_dir}")

        self.logger.info("=" * 60)
        self.logger.info("Method 1: BLAST + anicalc + aniclust")
        self.logger.info(f"  Input: {input_fasta}")
        self.logger.info(f"  Min ANI={self.min_ani}, Min tcov={self.min_tcov}, Min qcov={self.min_qcov}")
        self.logger.info("=" * 60)

        # ---------------------------
        # Step 1: 长度过滤
        # ---------------------------
        filtered_fasta = os.path.join(tmp_dir, f"{fasta_name}_len{self.min_length}.fasta")
        seq_lengths = _parse_fasta_lengths(input_fasta, self.min_length)
        self.logger.info(f"Step 1: Filtered {len(seq_lengths)} sequences >= {self.min_length}bp")

        if len(seq_lengths) < 2:
            self.logger.warning("  < 2 sequences after filtering, returning empty result")
            return ClusteringResult(method="blast_aniclust")

        _write_representatives(input_fasta, set(seq_lengths.keys()), filtered_fasta)

        # ---------------------------
        # Step 2: makeblastdb
        # ---------------------------
        db_path = os.path.join(tmp_dir, "temp_db")
        self.logger.info("Step 2: makeblastdb ...")
        try:
            subprocess.run([
                "makeblastdb", "-in", os.path.abspath(filtered_fasta),
                "-dbtype", "nucl", "-out", os.path.abspath(db_path)
            ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"  makeblastdb failed: {e.stderr.decode('utf-8', errors='replace') if e.stderr else 'unknown error'}")
            raise RuntimeError(f"makeblastdb failed. Check that BLAST+ is installed and working.") from e

        # ---------------------------
        # Step 3: blastn (全vs全)
        # ---------------------------
        blast_out = os.path.join(tmp_dir, "blast_output.tsv")
        self.logger.info("Step 3: blastn (all-vs-all) ...")
        try:
            subprocess.run([
                "blastn", "-query", os.path.abspath(filtered_fasta),
                "-db", os.path.abspath(db_path),
                "-outfmt", "6 std qlen slen",
                "-max_target_seqs", "10000",
                "-out", os.path.abspath(blast_out),
                "-num_threads", str(self.threads)
            ], check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"  blastn failed: {e}")
            raise RuntimeError(f"blastn failed.") from e

        if not os.path.getsize(blast_out):
            self.logger.warning("  BLAST produced no output, returning empty result")
            return ClusteringResult(method="blast_aniclust")

        # ---------------------------
        # Step 4: 计算配对ANI (内嵌anicalc)
        # ---------------------------
        ani_out = os.path.join(tmp_dir, "ani_output.tsv")
        self.logger.info("Step 4: Computing pairwise ANI (anicalc) ...")
        with open(blast_out) as fh_in, open(ani_out, 'w') as fh_out:
            fh_out.write("qname\ttname\tnum_alns\tpid\tqcov\ttcov\n")
            for alns in _yield_alignment_blocks(fh_in):
                alns = _prune_alns(alns)
                if not alns:
                    continue
                ani = _compute_ani(alns)
                qcov, tcov = _compute_cov(alns)
                fh_out.write(f"{alns[0]['qname']}\t{alns[0]['tname']}\t{len(alns)}\t{ani}\t{qcov}\t{tcov}\n")

        # ---------------------------
        # Step 5: 质心聚类 (内嵌aniclust)
        # ---------------------------
        self.logger.info("Step 5: Centroid clustering (aniclust) ...")
        clusters = _centroid_cluster(seq_lengths, ani_out,
                                     self.min_ani, self.min_tcov, self.min_qcov)

        # ---------------------------
        # Step 6: 提取代表序列
        # ---------------------------
        rep_fasta = os.path.join(output_dir, "representatives.fasta")
        clusters_tsv = os.path.join(output_dir, "clusters.tsv")
        reps = set(clusters.keys())

        self.logger.info("Step 6: Extracting representative sequences ...")
        _write_representatives(input_fasta, reps, rep_fasta)

        # 输出聚类表
        with open(clusters_tsv, 'w') as fh:
            fh.write("representative\tmembers\n")
            for rep, members in clusters.items():
                fh.write(f"{rep}\t{','.join(members)}\n")

        # 构建结果
        result = ClusteringResult(
            method="blast_aniclust",
            clusters=clusters,
            representatives=list(reps),
            n_clusters=len(clusters),
            n_singletons=sum(1 for v in clusters.values() if len(v) == 1),
            n_total_sequences=sum(len(v) for v in clusters.values()),
            cluster_sizes={k: len(v) for k, v in clusters.items()},
        )
        result.mean_cluster_size = (result.n_total_sequences / result.n_clusters
                                    if result.n_clusters > 0 else 0)

        self.logger.info(f"  Result: {result.n_clusters} clusters, "
                         f"{result.n_singletons} singletons, "
                         f"{result.n_total_sequences} total seqs")
        return result

    # ----------------------------------------------------------------
    # 方法2: dRep
    # ----------------------------------------------------------------

    def run_drep(self, input_path: str, output_dir: str) -> ClusteringResult:
        """
        运行 dRep 去冗余 (带病毒适应参数)

        参数来源: ViOTUcluster/Modules/drep_module.sh:36
        和 ViWrap/scripts/run_dRep.py:29

        Parameters
        ----------
        input_path : str
            输入的FASTA文件或基因组目录。若为目录则自动生成genome_list。
        output_dir : str
            dRep 输出目录

        Returns
        -------
        ClusteringResult
        """
        os.makedirs(output_dir, exist_ok=True)
        self.logger.info("=" * 60)
        self.logger.info("Method 2: dRep (virus-adapted parameters)")
        self.logger.info(f"  Input: {input_path}")
        self.logger.info(f"  sa={DREP_VIRUS_PARAMS['sa']}, nc={DREP_VIRUS_PARAMS['nc']}, "
                         f"min_length={self.min_length}")
        self.logger.info(f"  Score: sizeW=1, all other weights=0")
        self.logger.info("=" * 60)

        # 构建基因组列表
        if os.path.isfile(input_path):
            # 单FASTA: 先按序列拆分为单个基因组文件
            genome_dir = os.path.join(self.work_dir, "drep_genomes")
            os.makedirs(genome_dir, exist_ok=True)
            self._split_fasta_to_genomes(input_path, genome_dir)
            genome_list = os.path.join(self.work_dir, "genome_list.txt")
            genome_files = sorted(Path(genome_dir).glob("*.fasta"))
            with open(genome_list, 'w') as fh:
                for gf in genome_files:
                    fh.write(f"{gf.absolute()}\n")
            self.logger.info(f"  Split into {len(genome_files)} individual genomes")
        elif os.path.isdir(input_path):
            genome_dir = input_path
            genome_list = os.path.join(self.work_dir, "genome_list.txt")
            genome_files = sorted(Path(genome_dir).glob("*.fasta"))
            if not genome_files:
                genome_files = sorted(Path(genome_dir).glob("*.fna"))
            with open(genome_list, 'w') as fh:
                for gf in genome_files:
                    fh.write(f"{gf.absolute()}\n")
            self.logger.info(f"  Found {len(genome_files)} genomes in directory")
        else:
            raise FileNotFoundError(f"Input path not found: {input_path}")

        if not genome_files:
            self.logger.warning("  No genome files found, returning empty result")
            return ClusteringResult(method="drep")

        # 构建 dRep 命令 (病毒适应参数)
        drep_cmd = [
            "dRep", "dereplicate", output_dir,
            "-g", genome_list,
            "-p", str(self.threads),
            "-l", str(self.min_length),
            "-pa", str(DREP_VIRUS_PARAMS["pa"]),
            "-sa", str(DREP_VIRUS_PARAMS["sa"]),
            "-nc", str(DREP_VIRUS_PARAMS["nc"]),
            "-comW", str(DREP_VIRUS_PARAMS["comW"]),
            "-conW", str(DREP_VIRUS_PARAMS["conW"]),
            "-strW", str(DREP_VIRUS_PARAMS["strW"]),
            "-N50W", str(DREP_VIRUS_PARAMS["N50W"]),
            "-sizeW", str(DREP_VIRUS_PARAMS["sizeW"]),
            "-centW", str(DREP_VIRUS_PARAMS["centW"]),
        ]
        if DREP_VIRUS_PARAMS["ignore_genome_quality"]:
            drep_cmd.append("--ignoreGenomeQuality")
        if DREP_VIRUS_PARAMS["skip_plots"]:
            drep_cmd.append("--skip_plots")

        self.logger.info(f"  Command: {' '.join(drep_cmd)}")

        # 执行 dRep
        try:
            subprocess.run(drep_cmd, check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"  dRep failed: {e}")
            return ClusteringResult(method="drep")
        except FileNotFoundError:
            self.logger.error("  dRep not found in PATH. Please install: conda install -c bioconda drep")
            return ClusteringResult(method="drep")

        # 解析 dRep 输出
        self.logger.info("  Parsing dRep output ...")
        result = _parse_drep_clusters(output_dir, self.min_length)

        # 输出聚类TSV
        if result.clusters:
            clusters_tsv = os.path.join(output_dir, "clusters.tsv")
            with open(clusters_tsv, 'w') as fh:
                fh.write("representative\tmembers\n")
                for rep, members in result.clusters.items():
                    fh.write(f"{rep}\t{','.join(members)}\n")

        return result

    # ----------------------------------------------------------------
    # 方法3: 同时运行两种方法并对比
    # ----------------------------------------------------------------

    def run_both(self, input_fasta: str, output_dir: str) -> ComparisonReport:
        """
        同时运行 BLAST+aniclust 和 dRep, 输出对比报告

        Parameters
        ----------
        input_fasta : str
            输入FASTA文件
        output_dir : str
            总输出目录 (其下自动创建 blast_aniclust/ drep/ comparison/)

        Returns
        -------
        ComparisonReport
        """
        os.makedirs(output_dir, exist_ok=True)

        blast_outdir = os.path.join(output_dir, "blast_aniclust")
        drep_outdir = os.path.join(output_dir, "drep")
        comp_outdir = os.path.join(output_dir, "comparison")
        os.makedirs(comp_outdir, exist_ok=True)

        self.logger.info("=" * 70)
        self.logger.info("  Running BOTH methods for comparison")
        self.logger.info("=" * 70)

        report = ComparisonReport()

        # Run BLAST+aniclust
        t0 = time.time()
        report.blast_result = self.run_blast_aniclust(input_fasta, blast_outdir)
        blast_time = time.time() - t0

        # Run dRep
        t0 = time.time()
        report.drep_result = self.run_drep(input_fasta, drep_outdir)
        drep_time = time.time() - t0

        # ---- 对比分析 ----
        blast_reps = set(report.blast_result.representatives) if report.blast_result else set()
        drep_reps = set(report.drep_result.representatives) if report.drep_result else set()

        # 序列ID去后缀匹配 (dRep输出会加.fna/.fasta后缀)
        drep_reps_clean = set()
        for rid in drep_reps:
            drep_reps_clean.add(rid.rsplit('.', 1)[0])
        drep_reps = drep_reps_clean

        report.shared_representatives = blast_reps & drep_reps
        report.blast_only_representatives = blast_reps - drep_reps
        report.drep_only_representatives = drep_reps - blast_reps

        # 构建簇大小分布对比
        blast_sizes = (report.blast_result.cluster_sizes if report.blast_result else {})
        drep_sizes = (report.drep_result.cluster_sizes if report.drep_result else {})

        report.summary = {
            "blast_time_sec": round(blast_time, 1),
            "drep_time_sec": round(drep_time, 1),
            "blast_n_clusters": report.blast_result.n_clusters if report.blast_result else 0,
            "drep_n_clusters": report.drep_result.n_clusters if report.drep_result else 0,
            "blast_n_singletons": report.blast_result.n_singletons if report.blast_result else 0,
            "drep_n_singletons": report.drep_result.n_singletons if report.drep_result else 0,
            "n_shared_representatives": len(report.shared_representatives),
            "n_blast_only": len(report.blast_only_representatives),
            "n_drep_only": len(report.drep_only_representatives),
            "blast_mean_cluster_size": round(sum(blast_sizes.values()) / len(blast_sizes), 2) if blast_sizes else 0,
            "drep_mean_cluster_size": round(sum(drep_sizes.values()) / len(drep_sizes), 2) if drep_sizes else 0,
            "overlap_rate": round(len(report.shared_representatives) / max(len(blast_reps | drep_reps), 1) * 100, 1),
        }

        # 输出对比报告
        self._write_comparison_report(report, comp_outdir)

        return report

    # ----------------------------------------------------------------
    # 辅助方法
    # ----------------------------------------------------------------

    def _split_fasta_to_genomes(self, fasta_path: str, output_dir: str):
        """将单FASTA拆分为每个序列一个文件的基因组目录 (dRep需要)"""
        try:
            from Bio import SeqIO
            for rec in SeqIO.parse(fasta_path, "fasta"):
                if len(rec.seq) >= self.min_length:
                    safe_id = rec.id.replace('/', '_').replace('|', '_').replace(':', '_')
                    with open(os.path.join(output_dir, f"{safe_id}.fasta"), 'w') as fh:
                        SeqIO.write(rec, fh, "fasta")
        except ImportError:
            sid, seq = None, []
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith('>'):
                        if sid and len(''.join(seq)) >= self.min_length:
                            safe_id = sid.replace('/', '_').replace('|', '_').replace(':', '_')
                            with open(os.path.join(output_dir, f"{safe_id}.fasta"), 'w') as out:
                                out.write(f">{sid}\n{''.join(seq)}\n")
                        sid = line[1:].strip().split()[0]
                        seq = []
                    else:
                        seq.append(line.strip())
                if sid and len(''.join(seq)) >= self.min_length:
                    safe_id = sid.replace('/', '_').replace('|', '_').replace(':', '_')
                    with open(os.path.join(output_dir, f"{safe_id}.fasta"), 'w') as out:
                        out.write(f">{sid}\n{''.join(seq)}\n")

    def _write_comparison_report(self, report: ComparisonReport, output_dir: str):
        """输出对比报告文件"""
        # JSON summary
        summary_path = os.path.join(output_dir, "comparison_summary.json")
        with open(summary_path, 'w') as fh:
            json.dump(report.to_dict(), fh, indent=2, default=str)
        self.logger.info(f"  Summary: {summary_path}")

        # TSV对比表
        tsv_path = os.path.join(output_dir, "comparison_report.tsv")
        with open(tsv_path, 'w') as fh:
            fh.write("metric\tblast_aniclust\tdrep\n")
            s = report.summary
            fh.write(f"n_clusters\t{s['blast_n_clusters']}\t{s['drep_n_clusters']}\n")
            fh.write(f"n_singletons\t{s['blast_n_singletons']}\t{s['drep_n_singletons']}\n")
            fh.write(f"mean_cluster_size\t{s['blast_mean_cluster_size']}\t{s['drep_mean_cluster_size']}\n")
            fh.write(f"runtime_sec\t{s['blast_time_sec']}\t{s['drep_time_sec']}\n")
        self.logger.info(f"  TSV report: {tsv_path}")

        # Venn 文本摘要
        venn_path = os.path.join(output_dir, "venn_summary.txt")
        with open(venn_path, 'w') as fh:
            fh.write(f"BLAST+aniclust representatives: {s['blast_n_clusters']}\n")
            fh.write(f"dRep representatives:           {s['drep_n_clusters']}\n")
            fh.write(f"Shared:                          {len(report.shared_representatives)}\n")
            fh.write(f"BLAST-only:                      {len(report.blast_only_representatives)}\n")
            fh.write(f"dRep-only:                       {len(report.drep_only_representatives)}\n")
            fh.write(f"Overlap rate:                    {s['overlap_rate']}%\n")
        self.logger.info(f"  Venn summary: {venn_path}")

        # 打印终端摘要
        self.logger.info("")
        self.logger.info("=" * 50)
        self.logger.info("  COMPARISON SUMMARY")
        self.logger.info("=" * 50)
        self.logger.info(f"  BLAST+aniclust: {s['blast_n_clusters']} clusters "
                         f"({s['blast_n_singletons']} singletons) "
                         f"[{s['blast_time_sec']}s]")
        self.logger.info(f"  dRep:           {s['drep_n_clusters']} clusters "
                         f"({s['drep_n_singletons']} singletons) "
                         f"[{s['drep_time_sec']}s]")
        self.logger.info(f"  Shared reps:    {len(report.shared_representatives)}")
        self.logger.info(f"  Overlap rate:   {s['overlap_rate']}%")


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="病毒基因组去冗余评估 — 集成 BLAST+aniclust 和 dRep 两种方法",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 仅运行 BLAST+aniclust
  python dereplication.py -i contigs.fasta -m blast -o out_blast

  # 仅运行 dRep (需安装 dRep)
  python dereplication.py -i contigs.fasta -m drep -o out_drep

  # 同时运行两种方法并对比
  python dereplication.py -i contigs.fasta -m both -o out_comparison -t 16

dRep 参数说明 (病毒适应):
  所有参数从 ViOTUcluster 和 ViWrap 项目中提取, 专为病毒基因组优化:
  - --ignoreGenomeQuality   : 病毒基因组不完整, 不用CheckM
  - sizeW=1, 其他权重全为0   : 仅按基因组大小选代表序列
  - MASH预筛选 pa=0.8        : 80% ANI 粗聚类
  - ANIm精细聚类 sa=0.95     : 95% ANI 种水平阈值
  - nc=0.85                  : 85% 最小覆盖度
        """
    )
    parser.add_argument("-i", "--input", required=True,
                        help="输入 FASTA 文件或基因组目录")
    parser.add_argument("-o", "--output", default="./dereplication_out",
                        help="输出目录 (default: ./dereplication_out)")
    parser.add_argument("-m", "--method", choices=["blast", "drep", "both"],
                        default="both", help="去冗余方法 (default: both)")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="线程数 (default: 8)")
    parser.add_argument("--min_ani", type=float, default=95.0,
                        help="最小 ANI 阈值 (default: 95)")
    parser.add_argument("--min_tcov", type=float, default=85.0,
                        help="最小 target 覆盖率 (default: 85)")
    parser.add_argument("--min_qcov", type=float, default=0.0,
                        help="最小 query 覆盖率 (default: 0)")
    parser.add_argument("--min_length", type=int, default=3000,
                        help="最小序列/基因组长度 (default: 3000)")
    parser.add_argument("--keep_tmp", action="store_true",
                        help="保留临时文件")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="详细输出")

    args = parser.parse_args()

    # 日志
    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_file = os.path.join(args.output, "dereplication.log")
    os.makedirs(args.output, exist_ok=True)

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ]
    )
    logger = logging.getLogger(__name__)
    logger.info(f"dereplication.py started with method={args.method}")

    # 初始化
    work_dir = os.path.join(args.output, "_tmp") if not args.keep_tmp else os.path.join(args.output, "tmp")
    derep = Dereplicator(
        work_dir=work_dir,
        threads=args.threads,
        min_length=args.min_length,
        min_ani=args.min_ani,
        min_tcov=args.min_tcov,
        min_qcov=args.min_qcov,
    )

    # 执行
    if args.method == "blast":
        result = derep.run_blast_aniclust(args.input, args.output)
        logger.info(f"Done. {result.n_clusters} clusters found.")
    elif args.method == "drep":
        result = derep.run_drep(args.input, args.output)
        logger.info(f"Done. {result.n_clusters} clusters found.")
    else:  # both
        report = derep.run_both(args.input, args.output)
        s = report.summary
        logger.info(f"Comparison done. Overlap rate: {s['overlap_rate']}% "
                     f"({s['n_shared_representatives']} shared / "
                     f"{s['blast_n_clusters']}+{s['drep_n_clusters']} total)")

    if not args.keep_tmp:
        tmp_dir = os.path.join(args.output, "_tmp")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
            logger.info("Cleaned up temporary files.")


if __name__ == "__main__":
    main()
