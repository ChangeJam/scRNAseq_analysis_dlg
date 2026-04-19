#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize UMAP embeddings colored by cluster identity and QC metrics.

Produces 3 plots to assess clustering quality and QC effectiveness.
"""

import scanpy as sc
from pathlib import Path
import sys
import matplotlib.pyplot as plt

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

def main(config: dict) -> None:
    """Generate and save UMAP visualizations"""

    input_path = Path(config["paths"]["processed_dir"]) / "umap_result.h5ad"
    figdir = Path(config["paths"]["figures_dir"])
    fig1_path = figdir / "umap_leiden.png"
    fig2_path = figdir / "umap_ngenes.png"
    fig3_path = figdir / "umap_mt.png"

    figdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading UMAP-transformed data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Plot 1: UMAP colored by Leiden cluster
    print("Plotting UMAP by Leiden cluster...")
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='leiden', ax=ax1, show=False,
               title='Leiden Clusters')
    fig1.savefig(fig1_path, dpi=150, bbox_inches='tight')
    plt.close(fig1)
    print(f"    Saved: {fig1_path}")

    # Plot 2: UMAP colored by n_genes
    print("Plotting UMAP by n_genes...")
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='n_genes_by_counts', ax=ax2, show=False,
               title='Detected Genes per cell', cmap='viridis')
    fig2.savefig(fig2_path, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"    Saved: {fig2_path}")

    # Plot 3: UMAP colored by pct_counts_mt
    print("Plotting UMAP by pct_counts_mt...")
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='pct_counts_mt', ax=ax3, show=False,
               title='Mitochondrial Gene Percentage', cmap='RdYlBu_r')
    fig3.savefig(fig3_path, dpi=150, bbox_inches='tight')
    plt.close(fig3)
    print(f"    Saved: {fig3_path}")

    print("\nAll UMAP plots generated successfully!")

if __name__ == "__main__":
    config = get_config()
    main(config)
    print("UMAP visualization completed successfully.")