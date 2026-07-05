#!/usr/bin/env python3
"""
聚类去重雷达图 — 多方法 × 多指标对比
Dimensions: ARI, AMI, NMI, Homogeneity, Completeness, V-measure, Purity

用法:
  python plots/plot_clustering_radar.py \
      --input dedup_overall.tsv \
      --resource dedup_resource.tsv \
      --out_dir ./plots/
"""

import argparse
import os
import sys
from math import pi

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
           "#7E6148", "#8B5CF6", "#EC4899", "#F59E0B", "#06B6D4"]

CLUSTER_METRICS = ["ARI", "AMI", "NMI", "Homogeneity", "Completeness",
                   "V_measure", "Purity"]


def load_data(overall_path, resource_path=None):
    df = pd.read_csv(overall_path, sep="\t")
    if "Tool" not in df.columns:
        print(f"[ERROR] Missing 'Tool' column in {overall_path}")
        sys.exit(1)
    res = None
    if resource_path and os.path.exists(resource_path):
        res = pd.read_csv(resource_path, sep="\t")
    return df, res


def _norm(values, bigger_better=True):
    arr = np.array(values, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0)
    if arr.max() == arr.min():
        return np.ones_like(arr) * 0.5
    scaled = (arr - arr.min()) / (arr.max() - arr.min())
    if not bigger_better:
        scaled = 1.0 - scaled
    return scaled


def plot_radar(df_overall, df_res, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    tools = df_overall["Tool"].unique().tolist()
    if len(tools) < 2:
        print("[radar] Need at least 2 tools for comparison")
        return

    # 选择可用指标
    available = [m for m in CLUSTER_METRICS if m in df_overall.columns]
    if len(available) < 3:
        print("[radar] Need at least 3 metrics")
        return

    n = len(available)
    angles = [2 * pi * i / n for i in range(n)]
    angles += angles[:1]

    # 按 Tool 聚合均值
    agg = df_overall.groupby("Tool")[available].mean().reset_index()

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    for i, tool in enumerate(tools):
        row = agg[agg["Tool"] == tool]
        if row.empty:
            continue
        values = row[available].values[0].tolist()
        normed = _norm(values, bigger_better=True)
        normed = normed.tolist() + [normed[0]]

        color = PALETTE[i % len(PALETTE)]
        ax.fill(angles, normed, alpha=0.1, color=color)
        ax.plot(angles, normed, "o-", linewidth=2, color=color,
                label=tool, markersize=6)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(available, fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"],
                        fontsize=8, color="grey")
    ax.set_title("Clustering Dedup: Multi-Metric Radar",
                 fontsize=14, fontweight="bold", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=10)

    outpath = os.path.join(out_dir, "Fig_Clustering_Radar.png")
    plt.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[radar] {outpath}")


def main():
    p = argparse.ArgumentParser(description="聚类去重雷达图")
    p.add_argument("--input", required=True, help="dedup_overall.tsv")
    p.add_argument("--resource", default=None, help="dedup_resource.tsv (可选)")
    p.add_argument("--out_dir", default="plots/clustering")
    args = p.parse_args()

    if not os.path.exists(args.input):
        print(f"[ERROR] Input not found: {args.input}")
        sys.exit(1)

    df, res = load_data(args.input, args.resource)
    plot_radar(df, res, args.out_dir)


if __name__ == "__main__":
    main()
