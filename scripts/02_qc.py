#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quality control: filter cells using MAD-based outlier detection and Scrublet doublet removal.

Reads raw data, calculates QC metrics, removes outlier cells,
and saves filtered data.
"""

import pandas as pd
import scanpy as sc
import numpy as np
import scrublet as scr
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

def mad_outlier(adata, metric, n_thresh=5):
    """
    Identify outlier cells using Median Absolute Deviation.

    Args:
        adata: AnnData object
        metric: Column name in adata.obs to evaluate
        n_thresh: Number of MADs from median to define outlier boundary

    Returns:
        Boolean array where True marks an outlier cell
    """
    data = adata.obs[metric]
    #Calculate the median
    med = np.median(data)
    #Calculate the Median Absolute Deviation (MAD)
    mad = np.median(np.abs(data - med)) * 1.4826

    #Define the upper and lower bounds
    upper_bound = med + n_thresh * mad
    lower_bound = med - n_thresh * mad

    #Return a Boolean array marking cells that exceed the bounds as True
    return (data > upper_bound) | (data < lower_bound)

def detect_doublet(
        adata,
        expected_doublet_rate: float=0.08,
        random_state: int=0
) -> np.ndarray:
    """Run doublet detection on raw count data."""

    # ---初始化Scrublet---
    scrub = scr.Scrublet(
        adata.X,
        expected_doublet_rate=expected_doublet_rate,
        random_state=random_state
    )

    # ---运行双胞检测---
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # ---结果存入adata.obs---
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = pd.Categorical(predicted_doublets)


    # ---打印摘要---
    print(f"[Scrublet] expected doublet rate: {expected_doublet_rate:.2%}")
    print(f"[Scrublet] detected doublet rate: {np.mean(predicted_doublets):.2%}")
    print(f"[Scrublet] doublets predicted: {np.sum(predicted_doublets)}")

    # ---保存score的violin图---
    sc.pl.violin(
        adata,
        keys=["doublet_score"],
        groupby="predicted_doublet",
        jitter=0.4,
        save="_doublet_score.png"
    )

    return predicted_doublets
    

def main(config: dict) -> None:
    """Run QC filtering pipeline"""

    processed_dir = Path(config["paths"]["processed_dir"])
    figures_dir = Path(config["paths"]["figures_dir"])

    # ---参数---
    mt_genes = config["params"]["qc"]["mt_genes_complete"]
    mad_thresh = config["params"]["qc"]["mad_thresh"]
    expected_doublet_rate = config["params"]["qc"]["expected_doublet_rate"]
    scrublet_random_state = config["params"]["qc"]["scrublet_random_state"]

    # ---读取数据---
    input_path = processed_dir / "raw_data.h5ad"
    print(f"Reading data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # ---标记线粒体基因并计算QC指标---
    adata.var["mt"] = adata.var_names.isin(mt_genes)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # ---画QC小提琴图和散点图（过滤前）---
    qc_plot_dir = figures_dir / "qc_plots"
    qc_plot_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(qc_plot_dir)

    sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_before_filter.png'
    )

    sc.pl.scatter(
    adata, 
    x='total_counts', 
    y='n_genes_by_counts', 
    color='pct_counts_mt',
    save='_before_filter.png'
    )

    # ---Scrublet双胞检测---
    print("\n---Running Scrublet doublet detection---")
    predicted_doublets = detect_doublet(
        adata,
        expected_doublet_rate=expected_doublet_rate,
        random_state=scrublet_random_state
    )

    # ---MAD离群值过滤---
    print(f"n_obs(before filter): {adata.n_obs}")

    outlier_genes = mad_outlier(adata, "n_genes_by_counts", n_thresh=mad_thresh)
    outlier_counts = mad_outlier(adata, "total_counts", n_thresh=mad_thresh)
    outlier_mt = mad_outlier(adata, "pct_counts_mt", n_thresh=mad_thresh)

    # ---合并MAD outlier和scrublet doublet
    outlier_mask = outlier_genes | outlier_mt | outlier_counts |predicted_doublets

    adata = adata[~outlier_mask, :].copy()

    print(f"n_obs(after filter): {adata.n_obs}")

    # ---画QC小提琴图和散点图（过滤后）---
    sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_after_filter.png'
    )

    sc.pl.scatter(
    adata, 
    x='total_counts', 
    y='n_genes_by_counts', 
    color='pct_counts_mt',
    save='_after_filter.png'
    )

    # ---保存过滤后的数据---
    output_path = processed_dir / "filtered_data.h5ad"
    adata.write_h5ad(output_path)
    print(f"QC done, filtered data saved to: {output_path}")


if __name__ == "__main__":
    config = get_config()
    main(config)
    print("Step 02 completed: QC filtered done (MAD + Scrublet).")