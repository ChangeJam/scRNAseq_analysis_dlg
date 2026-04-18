#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quality control: filter cells using MAD-based outlier detection.

Reads raw data, calculates QC metrics, removes outlier cells,
and saves filtered data.
"""

import scanpy as sc
import numpy as np
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


def main(config: dict) -> None:
    """Run QC filtering pipeline"""

    processed_dir = Path(config["paths"]["processed_dir"])
    figures_dir = Path(config["paths"]["figures_dir"])

    # ---参数---
    mt_genes = config["params"]["qc"]["mt_genes"]
    mad_thresh = config["params"]["qc"]["mad_thresh"]

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
    save='_before_filter.pdf'
    )

    sc.pl.scatter(
    adata, 
    x='total_counts', 
    y='n_genes_by_counts', 
    color='pct_counts_mt',
    save='_before_filter.pdf'
    )

    # ---MAD离群值过滤---
    print(f"n_obs(before filter): {adata.n_obs}")

    outlier_genes = mad_outlier(adata, "n_genes_by_counts", n_thresh=mad_thresh)
    outlier_counts = mad_outlier(adata, "total_counts", n_thresh=mad_thresh)
    outlier_mt = mad_outlier(adata, "pct_counts_mt", n_thresh=mad_thresh)

    adata = adata[~(outlier_genes | outlier_counts | outlier_mt), :].copy()

    print(f"n_obs(after filter): {adata.n_obs}")

    # ---保存过滤后的数据---
    output_path = processed_dir / "filtered_data.h5ad"
    adata.write_h5ad(output_path)
    print(f"QC done, filtered data saved to: {output_path}")



if __name__ == "__main__":
    config = get_config()
    main(config)
    print("Step 02 completed: QC filtered done.")