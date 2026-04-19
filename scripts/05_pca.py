#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA Dimension Reduction.

Performs PCA on HVG-subset data and saves the result.
"""
import scanpy as sc
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

def main(config: dict) -> None:
    """
    Execute PCA and save the resulting AnnData object.
    """
    input_path = Path(config["paths"]["processed_dir"]) / "hvg_subset.h5ad"
    output_path = Path(config["paths"]["processed_dir"]) / "pca_result.h5ad"
    n_comps = config["params"]["dim_reduction"]["n_pcs"]

    print(f"Loading HVG-subset data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    print(f"Running PCA to compute {n_comps} principal components...")
    sc.pp.pca(adata, n_comps=n_comps)

    # ---检查输出矩阵形状---
    print(f"PCA output matrix shape (Cells x PCs): {adata.obsm['X_pca'].shape}")
    print(f"PC loadings matrix shape (Genes x PCs): {adata.varm['PCs'].shape}")

    print(f"Saving PCA-transformed data to: {output_path}")
    adata.write_h5ad(output_path)


if __name__ == "__main__":
    config = get_config()
    main(config)
    print("PCA computation completed successfully.")