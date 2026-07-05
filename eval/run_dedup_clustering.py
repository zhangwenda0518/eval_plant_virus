#!/usr/bin/env python3
"""
序列去冗余聚类完整评估管道
借鉴 ViralClust 思路，集成 5 种聚类工具 + 金标准 ARI/AMI 评估
工具: CD-HIT, MMseqs2, VCLUST (VSEARCH), SUMACLUST, HDBSCAN (UMAP+密度)

用法:
  python run_dedup_clustering.py \
      --input step2_dedup_fragments.fasta \
      --outdir step3_dedup_cluster/ \
      --threads 30 --seed 42
"""

import argparse, os, sys, subprocess, shutil, glob, re, gzip, time
from collections import defaultdict
import numpy as np
import pandas as pd

try:
    import resource
    HAS_RESOURCE = True
except ImportError:
    HAS_RESOURCE = False
from sklearn.metrics import (adjusted_rand_score, adjusted_mutual_info_score,
                              homogeneity_score, completeness_score,
                              v_measure_score, normalized_mutual_info_score)


# ══════════════════════════════════════
# Part 0: 工具检查
# ══════════════════════════════════════

# 默认搜索路径 + 常见软件安装位置
_EXTRA_PATHS = [
    os.path.expanduser("~/biosoft/binary"),
    os.path.expanduser("~/sysoft/Python-3.7.10/bin"),
    os.path.expanduser("~/mambaforge/bin"),
    os.path.expanduser("~/bin"),
]

def _find_exe(name, extra_paths=None):
    """查找可执行文件：先 PATH，再额外路径"""
    found = shutil.which(name)
    if found: return found
    for d in (extra_paths or []):
        p = os.path.join(d, name)
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    # 也搜索所有默认额外路径
    for d in _EXTRA_PATHS:
        if os.path.isdir(d):
            for f in os.listdir(d):
                fp = os.path.join(d, f)
                if os.path.isfile(fp) and os.access(fp, os.X_OK):
                    if f == name or f.startswith(name):
                        return fp
    return None


def check_tools(tool_paths=None):
    """检测可用工具 {名称: 可执行文件路径}"""
    available = {}
    for name, default_exe in [('CD-HIT', 'cd-hit-est'), ('MMseqs2', 'mmseqs'),
                               ('VCLUST', 'vclust'),  # 三步: prefilter+align+cluster ('SUMACLUST', 'sumaclust'),
                               ('HDBSCAN', 'python3')]:
        # 用户指定路径优先
        exe = (tool_paths or {}).get(name, default_exe)
        found = _find_exe(exe, [os.path.dirname(exe)] if '/' in exe else None)
        if found:
            available[name] = found
    return available


# ── 资源追踪 ──
_RESOURCE_LOG = []

def tracked_run(tool_name, cmd, **kw):
    """执行子进程并记录墙钟时间、峰值内存"""
    t0 = time.time()
    # 尝试用 /usr/bin/time -v 获取精确内存
    if HAS_RESOURCE:
        p = subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **kw)
        t1 = time.time()
        mem_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    else:
        p = subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **kw)
        t1 = time.time()
        mem_kb = 0

    wall_s = t1 - t0
    mem_mb = mem_kb / 1024.0 if mem_kb > 0 else 0
    mem_str = f", {mem_mb:.0f} MB" if mem_mb > 0 else ""
    print(f"    ⏱ {wall_s:.1f}s{mem_str}")
    _RESOURCE_LOG.append({'Tool': tool_name, 'Wall_Time_s': round(wall_s, 1),
                          'Max_RSS_MB': round(mem_mb, 1)})
    return p


def save_resource_log(out_dir):
    if _RESOURCE_LOG:
        df = pd.DataFrame(_RESOURCE_LOG)
        df.to_csv(os.path.join(out_dir, "dedup_resource.tsv"), sep='\t', index=False)
        print(f"   📊 dedup_resource.tsv")
        return df
    return pd.DataFrame()


def plot_resource(df, out_dir):
    if df.empty: return
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    # 时间
    sns.barplot(data=df, x='Tool', y='Wall_Time_s', hue='Tool', palette=MY_PALETTE, legend=False, ax=axes[0])
    axes[0].set_ylabel('Wall Time (s)'); axes[0].set_title('Runtime', fontweight='bold')
    for p in axes[0].patches:
        axes[0].annotate(f'{p.get_height():.0f}s', (p.get_x()+p.get_width()/2, p.get_height()),
                        ha='center', va='bottom', fontsize=8)
    if df['Max_RSS_MB'].max() > 0:
        sns.barplot(data=df, x='Tool', y='Max_RSS_MB', hue='Tool', palette=MY_PALETTE, legend=False, ax=axes[1])
        axes[1].set_ylabel('Max RSS (MB)'); axes[1].set_title('Peak Memory', fontweight='bold')
        for p in axes[1].patches:
            axes[1].annotate(f'{p.get_height():.0f}MB', (p.get_x()+p.get_width()/2, p.get_height()),
                            ha='center', va='bottom', fontsize=8)
    else:
        axes[1].set_visible(False)
    fig.suptitle('Resource Consumption by Clustering Tool', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "dedup_resource.png"), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"   📊 dedup_resource.png")


# ══════════════════════════════════════
# Part 1: 运行聚类工具 → 统一 .clstr 格式
# ══════════════════════════════════════

def run_cdhit(input_fa, out_dir, threads, min_id, wordsize=10, exe="cd-hit-est"):
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = os.path.join(out_dir, "cdhit")
    cmd = [exe, "-i", input_fa, "-o", out_prefix,
           "-c", str(min_id), "-n", str(wordsize), "-d", "0",
           "-M", "32000", "-T", str(threads)]
    tracked_run('cdhit', cmd)
    # CD-HIT 原生就是 .clstr → 直接返回
    clstr_file = out_prefix + ".clstr"
    return parse_cdhit_clstr(clstr_file) if os.path.exists(clstr_file) else {}


def run_mmseqs(input_fa, out_dir, threads, min_id, exe="mmseqs"):
    os.makedirs(out_dir, exist_ok=True)
    tmp = os.path.join(out_dir, "tmp")
    prefix = os.path.join(out_dir, "mmseqs")
    # linclust: 比 easy-cluster 更快，--cov-mode 1 (target coverage) 适合不等长片段
    cmd = [exe, "easy-linclust", input_fa, prefix, tmp,
           "--min-seq-id", str(min_id), "-c", "0.8", "--cov-mode", "1",
           "--threads", str(threads), "--kmer-per-seq-scale", "0.3"]
    tracked_run('mmseqs', cmd)

    # MMseqs2 输出 rep \t member → 转为 .clstr 格式
    tsv_file = prefix + "_cluster.tsv"
    clusters = defaultdict(set)
    for line in open(tsv_file):
        rep, member = line.strip().split('\t')
        clusters[rep].add(rep)
        clusters[rep].add(member)
    return dict(clusters)


def run_vclust(input_fa, out_dir, threads, min_id, exe="vclust", qcov=0.3):
    """VCLUST: prefilter → align → cluster (Leiden vOTU, §6.2)"""
    os.makedirs(out_dir, exist_ok=True)

    # 1. prefilter: k-mer 预过滤
    fltr = os.path.join(out_dir, "fltr.txt")
    cmd1 = [exe, "prefilter", "-i", input_fa, "-o", fltr,
            "--min-ident", str(min_id), "--min-kmers", "15"]
    tracked_run('vclust', cmd1)

    # 2. align: 成对 ANI, 输出预过滤减小文件
    ani = os.path.join(out_dir, "ani.tsv")
    cmd2 = [exe, "align", "-i", input_fa, "-o", ani, "--filter", fltr,
            "--filter-threshold", str(min_id), "--outfmt", "lite",
            "--out-ani", str(min_id), "--out-qcov", str(qcov)]
    tracked_run('vclust', cmd2)

    # 3. cluster: Leiden 社区发现 (短片段靠传递连通)
    ids = ani + ".ids.tsv"
    clusters_out = os.path.join(out_dir, "vclust_clusters.tsv")
    cmd3 = [exe, "cluster", "-i", ani, "-o", clusters_out, "--ids", ids,
            "--algorithm", "leiden", "--metric", "ani",
            "--ani", str(min_id), "--qcov", str(qcov)]
    try:
        tracked_run('vclust', cmd3)
    except Exception as e:
        print(f"    ⚠ VCLUST cluster failed: {e}")
        return {}

    # 解析: object\tcluster_id
    clusters = defaultdict(set)
    for line in open(clusters_out):
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            clusters[parts[1]].add(parts[0])
    return dict(clusters)


def run_sumaclust(input_fa, out_dir, threads, min_id, exe="sumaclust"):
    """SUMACLUST — 基于后缀数组的快速聚类 (类似 SWARM)"""
    os.makedirs(out_dir, exist_ok=True)
    otu_out = os.path.join(out_dir, "sumaclust.otu")
    diff = int(10 * (1 - min_id) * 100)
    cmd = [exe, "-i", input_fa, "-o", otu_out,
           "-t", str(diff), "-p", str(threads)]
    tracked_run('SUMACLUST', cmd)
    return parse_sumaclust_otu(otu_out)


def run_hdbscan(input_fa, out_dir, threads, seed, exe="python3"):
    """HDBSCAN — UMAP 降维 + 密度聚类 (非参数，适合高度多态病毒)"""
    os.makedirs(out_dir, exist_ok=True)
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "viralclust", "bin", "hdbscan_virus.py")
    # 如果本地没有，创建一个简化版
    if not os.path.exists(script):
        script = _write_hdbscan_script(out_dir)

    cmd = [exe, script, "-i", input_fa, "-o", out_dir,
           "-t", str(threads), "-s", str(seed)]
    tracked_run('HDBSCAN', cmd)
    clstr_file = os.path.join(out_dir, "hdbscan_clusters.tsv")
    return parse_hdbscan_tsv(clstr_file) if os.path.exists(clstr_file) else {}


def _write_hdbscan_script(out_dir):
    """生成 HDBSCAN 聚类脚本"""
    script_path = os.path.join(out_dir, "_hdbscan_cluster.py")
    with open(script_path, 'w') as f:
        f.write("""#!/usr/bin/env python3
import argparse, os, sys
import numpy as np
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
import umap
import hdbscan

def kmer_features(fasta, k=5):
    ids, feats = [], []
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper()
        kmers = [seq[i:i+k] for i in range(len(seq)-k+1) if 'N' not in seq[i:i+k]]
        ids.append(rec.id.split()[0])
        feats.append(' '.join(kmers))
    vec = TfidfVectorizer(max_features=5000, sublinear_tf=True)
    X = vec.fit_transform(feats)
    svd = TruncatedSVD(n_components=min(50, X.shape[1]-1), random_state=42)
    X_red = svd.fit_transform(X)
    u = umap.UMAP(n_components=5, random_state=42, n_jobs=1)
    emb = u.fit_transform(X_red)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1, core_dist_n_jobs=1)
    labels = clusterer.fit_predict(emb)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hdbscan_clusters.tsv")
    with open(out, 'w') as fo:
        for sid, lb in zip(ids, labels):
            fo.write(f"{sid}\\t{lb}\\n")
    print(f"HDBSCAN: {len(set(labels))} clusters, {list(labels).count(-1)} noise")
""")
    return script_path


# ══════════════════════════════════════
# Part 2: 格式解析
# ══════════════════════════════════════

def parse_cdhit_clstr(clstr_file):
    clusters = defaultdict(set)
    current = None
    for line in open(clstr_file):
        if line.startswith('>'):
            current = line.strip().split()[1]
        elif current and '>' in line:
            seq_id = line.split('>')[1].split('...')[0].strip()
            clusters[current].add(seq_id)
    return dict(clusters)


def parse_vclust_uc(uc_file):
    clusters = defaultdict(set)
    for line in open(uc_file):
        p = line.strip().split('\t')
        if p[0] in ('C', 'H'):
            clusters[p[1]].add(p[8])
    return dict(clusters)


def parse_sumaclust_otu(otu_file):
    """SUMACLUST OTU 格式: cluster_id\\tmember1\\tmember2\\t..."""
    clusters = defaultdict(set)
    for line in open(otu_file):
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            cid = parts[0]
            for m in parts[1:]:
                clusters[cid].add(m)
    return dict(clusters)


def parse_hdbscan_tsv(tsv_file):
    """HDBSCAN 输出: seq_id\\tcluster_label"""
    clusters = defaultdict(set)
    for line in open(tsv_file):
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            clusters[parts[1]].add(parts[0])
    return dict(clusters)


# ══════════════════════════════════════
# Part 3: 金标准解析（同 eval_dedup_clustering.py）
# ══════════════════════════════════════

def parse_gold(fasta_path):
    labels, attrs = {}, {}
    for line in open(fasta_path):
        if not line.startswith('>'): continue
        sid = line[1:].strip().split()[0]
        parts = sid.split('_')
        if parts[0] == 'HOST':
            labels[sid] = '_'.join(parts[:2])
            mut = 'host'; lf = '?'
        else:
            labels[sid] = parts[0]
            mut = '0%'; lf = '?'
        for p in parts:
            if p.startswith('mut'): mut = p.replace('mut','').replace('pct','%')
            if p.startswith('len'): lf = p.replace('len','').replace('pct','%')
        attrs[sid] = {'mutation': mut, 'length_fraction': lf}
    return labels, attrs


def compute_metrics(true_labels, pred_clusters):
    all_ids = list(true_labels.keys())
    y_true = [true_labels[i] for i in all_ids]
    pmap = {}
    for cname, members in pred_clusters.items():
        for m in members:
            pmap[m] = cname
    y_pred = [pmap.get(i, f'__u_{hash(i)}') for i in all_ids]
    return {
        'ARI': adjusted_rand_score(y_true, y_pred),
        'AMI': adjusted_mutual_info_score(y_true, y_pred),
        'NMI': normalized_mutual_info_score(y_true, y_pred),
        'Homogeneity': homogeneity_score(y_true, y_pred),
        'Completeness': completeness_score(y_true, y_pred),
        'V_measure': v_measure_score(y_true, y_pred),
        'n_seqs': len(all_ids), 'n_species': len(set(y_true)),
        'n_clusters': len(pred_clusters),
    }


def stratified_metrics(true_labels, attrs, clusters, stratify_by):
    results = []
    strata = set(a[stratify_by] for a in attrs.values())
    for val in sorted(strata):
        sids = {i for i, a in attrs.items() if a[stratify_by] == val}
        sub_true = {i: true_labels[i] for i in sids if i in true_labels}
        sub_pred = {}; sids_set = set(sub_true.keys())
        for cn, ms in clusters.items():
            fm = ms & sids_set
            if fm: sub_pred[cn] = fm
        if sub_true and sub_pred:
            m = compute_metrics(sub_true, sub_pred)
            m[stratify_by] = val; results.append(m)
    return pd.DataFrame(results)


# ══════════════════════════════════════
# Part 4: 可视化
# ══════════════════════════════════════

MY_PALETTE = {'CD-HIT':'#4C72B0','MMseqs2':'#55A868','VCLUST':'#DD8452',
              'SUMACLUST':'#8B5CF6','HDBSCAN':'#E69F00'}
MUT_ORDER = ['0%','5%','host']

def plot_results(overall, mut_all, len_all, out_dir, min_id):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(style='ticks', font_scale=1.05)

    for metric, yl in [('ARI','Adjusted Rand Index'), ('V_measure','V-measure')]:
        order = overall.sort_values(metric, ascending=False)['Tool']
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.barplot(data=overall, x='Tool', y=metric, order=order, hue='Tool',
                    palette=MY_PALETTE, legend=False, ax=ax)
        ax.set_ylim(0, 1.02); ax.set_ylabel(yl)
        ax.set_title(f"{yl} (id={min_id})", fontweight='bold')
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"bar_{metric}.png"), dpi=200, bbox_inches='tight')
        plt.close()

        if not mut_all.empty:
            fig, ax = plt.subplots(figsize=(8, 5))
            sns.boxplot(data=mut_all, x='mutation', y=metric, hue='Tool',
                        order=MUT_ORDER, hue_order=[t for t in overall['Tool'] if t in mut_all['Tool'].unique()],
                        palette=MY_PALETTE, ax=ax)
            ax.set_ylim(0, 1.02); ax.set_ylabel(yl); ax.set_title(f"{yl} by Mutation", fontweight='bold')
            ax.legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"box_{metric}_by_mutation.png"), dpi=200, bbox_inches='tight')
            plt.close()

        if not len_all.empty:
            fig, ax = plt.subplots(figsize=(10, 5))
            sns.lineplot(data=len_all, x='length_fraction', y=metric, hue='Tool',
                        hue_order=[t for t in overall['Tool'] if t in len_all['Tool'].unique()],
                        palette=MY_PALETTE, markers=True, ax=ax)
            ax.set_ylim(0, 1.02); ax.set_ylabel(yl)
            ax.set_title(f"{yl} by Fragment Length", fontweight='bold'); ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"line_{metric}_by_length.png"), dpi=200, bbox_inches='tight')
            plt.close()


# ══════════════════════════════════════
# Main
# ══════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="去冗余聚类完整评估管道")
    parser.add_argument("--input", required=True, help="模拟片段 FASTA")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--min-id", type=float, default=0.90)
    parser.add_argument("--qcov", type=float, default=0.30,
                        help="覆盖度阈值 (query-aligned fraction, 0-1, 片段数据建议 0.3-0.5)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-run", action="store_true", help="跳过聚类运行，仅评估")
    parser.add_argument("--mutation", default=None,
                        help="仅评估指定突变率 (0 或 5)，如 --mutation 0")
    parser.add_argument("--tools", default="cdhit,mmseqs,vclust,sumaclust,hdbscan",
                        help="逗号分隔 (default: all 5)")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    tool_list = [t.strip() for t in args.tools.split(',')]

    # 如果指定了突变率过滤，先拆分 FASTA
    input_fa = args.input
    if args.mutation is not None:
        mut_tag = f"mut{args.mutation}pct"
        input_fa = os.path.join(args.outdir, f"_filtered_{mut_tag}.fasta")
        if not os.path.exists(input_fa):
            print(f"  Filtering {mut_tag} fragments...")
            n_in, n_out = 0, 0
            with open(args.input) as fin, open(input_fa, 'w') as fout:
                write = False
                for line in fin:
                    if line.startswith('>'):
                        write = mut_tag in line
                    if write:
                        fout.write(line)
                        n_out += 1
                    n_in += 1
            print(f"    → {n_out} lines written")
        else:
            print(f"  Using cached filter: {input_fa}")

    print("═" * 50)
    print("  去冗余聚类评估管道")
    print("═" * 50)
    print(f"  Input:  {input_fa} ({grep_count(input_fa)} sequences)")
    print(f"  Min ID: {args.min_id}")
    if args.mutation:
        print(f"  Mutation: {args.mutation}%")
    print(f"  Tools:  {', '.join(tool_list)}")
    print("═" * 50)

    # ── 运行聚类 ──
def parse_aniclust_clusters(tsv_file):
    """BLAST+aniclust / dRep clusters.tsv: rep\\tmember1,member2,..."""
    clusters = defaultdict(set)
    for line in open(tsv_file):
        if line.startswith('representative'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            rep = parts[0]
            members = [m.strip() for m in parts[1].split(',') if m.strip()]
            clusters[rep].add(rep)
            for m in members:
                clusters[rep].add(m)
    return dict(clusters)


def run_blast_aniclust(input_fa, out_dir, threads, min_id, qcov=0.3):
    """BLAST+aniclust: all-vs-all BLASTn → weighted ANI → centroid clustering"""
    os.makedirs(out_dir, exist_ok=True)
    derep_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "metrics", "dereplication.py")
    if not os.path.exists(derep_script):
        derep_script = "dereplication.py"  # assume in PATH

    ani_pct = int(min_id * 100)
    tcov = int(qcov * 100)
    cmd = [sys.executable, derep_script, "-i", input_fa, "-o", out_dir,
           "-m", "blast", "-t", str(threads),
           "--min_ani", str(ani_pct), "--min_tcov", str(tcov),
           "--min_length", "200"]
    tracked_run('blast', cmd)
    tsv = os.path.join(out_dir, "clusters.tsv")
    return parse_aniclust_clusters(tsv) if os.path.exists(tsv) else {}


def run_drep_cluster(input_fa, out_dir, threads, min_id, qcov=0.3):
    """dRep: MASH prefilter + ANIm clustering (病毒适应参数)"""
    os.makedirs(out_dir, exist_ok=True)
    derep_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "metrics", "dereplication.py")
    if not os.path.exists(derep_script):
        derep_script = "dereplication.py"

    ani_pct = int(min_id * 100)
    tcov = int(qcov * 100)
    cmd = [sys.executable, derep_script, "-i", input_fa, "-o", out_dir,
           "-m", "drep", "-t", str(threads),
           "--min_ani", str(ani_pct), "--min_tcov", str(tcov),
           "--min_length", "200"]
    tracked_run('drep', cmd)
    tsv = os.path.join(out_dir, "clusters.tsv")
    return parse_aniclust_clusters(tsv) if os.path.exists(tsv) else {}


RUNNERS = {
        'cdhit':    lambda: run_cdhit(input_fa, os.path.join(args.outdir, 'cdhit'), args.threads, args.min_id),
        'mmseqs':   lambda: run_mmseqs(input_fa, os.path.join(args.outdir, 'mmseqs'), args.threads, args.min_id),
        'vclust':   lambda: run_vclust(input_fa, os.path.join(args.outdir, 'vclust'), args.threads, args.min_id, qcov=args.qcov),
        'sumaclust':lambda: run_sumaclust(input_fa, os.path.join(args.outdir, 'sumaclust'), args.threads, args.min_id),
        'hdbscan':  lambda: run_hdbscan(input_fa, os.path.join(args.outdir, 'hdbscan'), args.threads, args.seed),
        'blast':    lambda: run_blast_aniclust(input_fa, os.path.join(args.outdir, 'blast'), args.threads, args.min_id, qcov=args.qcov),
        'drep':     lambda: run_drep_cluster(input_fa, os.path.join(args.outdir, 'drep'), args.threads, args.min_id, qcov=args.qcov),
    }

    tools = {}
    for tool in tool_list:
        if args.skip_run:
            print(f"  [SKIP] {tool} (--skip-run)")
            continue
        if tool not in RUNNERS:
            print(f"  [SKIP] {tool} not supported")
            continue
        print(f"  [{tool}] clustering...")
        try:
            tools[tool] = RUNNERS[tool]()
            print(f"    → {len(tools[tool])} clusters")
        except Exception as e:
            print(f"    ✗ FAILED: {e}")

    if not tools:
        print("[ERROR] No clustering results"); return

    # ── 金标准解析 ──
    true_labels, attrs = parse_gold(input_fa)

    # ── 评估 ──
    overall_rows, mut_rows, len_rows = [], [], []
    for name, clusters in tools.items():
        m = compute_metrics(true_labels, clusters)
        m['Tool'] = name; overall_rows.append(m)
        mdf = stratified_metrics(true_labels, attrs, clusters, 'mutation')
        if not mdf.empty: mdf['Tool'] = name; mut_rows.append(mdf)
        ldf = stratified_metrics(true_labels, attrs, clusters, 'length_fraction')
        if not ldf.empty: ldf['Tool'] = name; len_rows.append(ldf)

    overall = pd.DataFrame(overall_rows)
    print(f"\n  🏆 Overall (id={args.min_id}):")
    for _, r in overall.iterrows():
        print(f"    {r['Tool']:12s}  ARI={r['ARI']:.4f}  AMI={r['AMI']:.4f}  "
              f"V={r['V_measure']:.4f}  NMI={r['NMI']:.4f}  n_clust={r['n_clusters']}")

    overall.to_csv(os.path.join(args.outdir, "dedup_overall.tsv"), sep='\t', index=False)
    mut_all = pd.concat(mut_rows) if mut_rows else pd.DataFrame()
    len_all = pd.concat(len_rows) if len_rows else pd.DataFrame()
    if not mut_all.empty:
        mut_all.to_csv(os.path.join(args.outdir, "dedup_by_mutation.tsv"), sep='\t', index=False)
    if not len_all.empty:
        len_all.to_csv(os.path.join(args.outdir, "dedup_by_length.tsv"), sep='\t', index=False)

    plot_results(overall, mut_all, len_all, args.outdir, args.min_id)
    # 资源评估
    res_df = save_resource_log(args.outdir)
    if not res_df.empty:
        plot_resource(res_df, args.outdir)
    print(f"\n✅ Done. Output: {args.outdir}")


def grep_count(fasta):
    n = 0
    for line in open(fasta):
        if line.startswith('>'): n += 1
    return n


if __name__ == "__main__":
    main()
