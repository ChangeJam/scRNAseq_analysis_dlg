#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualize UMAP embeddings colored by cluster identity and QC metrics.
Produces three plots to assess clustering quality and QC effectiveness.
"""

import scanpy as sc
import os
import matplotlib.pyplot as plt

# Output directory for figures
FIG_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/figures"

def plot_umap(input_path: str) -> None:
    """Generate and save UMAP visualizations"""

    os.makedirs(FIG_DIR, exist_ok=True)

    print(f"Loading UMAP-transformed data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Plot 1: UMAP colored by Leiden cluster
    print("Plotting UMAP by Leiden cluster...")
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='leiden', ax=ax1, show=False,
               title='Leiden Clusters')
    fig1.savefig(os.path.join(FIG_DIR, "umap_leiden.png"), dpi=150, bbox_inches='tight')
    plt.close(fig1)
    print(f"    Saved: {os.path.join(FIG_DIR, 'umap_leiden.png')}")

    # Plot 2: UMAP colored by n_genes
    print("Plotting UMAP by n_genes...")
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='n_genes_by_counts', ax=ax2, show=False,
               title='Detected Genes per cell', cmap='viridis')
    fig2.savefig(os.path.join(FIG_DIR, "umap_ngenes.png"), dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"    Saved: {os.path.join(FIG_DIR, 'umap_ngenes.png')}")

    # Plot 3: UMAP colored by pct_counts_mt
    print("Plotting UMAP by pct_counts_mt...")
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='pct_counts_mt', ax=ax3, show=False,
               title='Mitochondrial Gene Percentage', cmap='RdYlBu_r')
    fig3.savefig(os.path.join(FIG_DIR, "umap_mt.png"), dpi=150, bbox_inches='tight')
    plt.close(fig3)
    print(f"    Saved: {os.path.join(FIG_DIR, 'umap_mt.png')}")

    print("\nAll UMAP plots generated successfully!")

if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "umap_result.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)

    plot_umap(input_path)