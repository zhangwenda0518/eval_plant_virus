#!/usr/bin/env python3
"""统一绘图工具模块 — 颜色方案、字体、保存格式"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 颜色方案
PALETTE = {
    "Salmon": "#4C72B0",
    "Bowtie2": "#DD8452",
    "Minimap2": "#55A868",
    "Kraken2": "#C44E52",
    "Diamond": "#8B5CF6",
    "MEGAHIT": "#4C72B0",
    "rnaviralspades": "#DD8452",
    "Penguin": "#55A868",
    "RefineC_Merge": "#C44E52",
    "V1_no_verify": "#C44E52",
    "V2_uniref90": "#DD8452",
    "V3_adversarial": "#4C72B0",
    "MMseqs2": "#4C72B0",
    "VITAP": "#DD8452",
    "ACVirus": "#55A868",
    "Integrated": "#C44E52",
    "D0": "#95A5A6",
    "D1": "#4C72B0",
    "D2": "#DD8452",
    "D3": "#55A868",
    "D4": "#C44E52",
    "positive": "#2ECC71",
    "negative_A": "#95A5A6",
    "negative_B": "#E74C3C",
    "negative_C": "#F39C12",
}

# 样式
def set_style():
    sns.set_theme(style="ticks", font_scale=1.2, context="paper")
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'axes.titlesize': 14,
        'axes.labelsize': 13,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 11,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })

def get_color(name):
    return PALETTE.get(name, "#333333")

def save_fig(fig, path, formats=('png',)):
    for fmt in formats:
        fig.savefig(f"{path}.{fmt}", dpi=300, bbox_inches='tight')
    plt.close(fig)

def save_default(fig, path):
    save_fig(fig, path)
    print(f"  Saved: {path}.png")
