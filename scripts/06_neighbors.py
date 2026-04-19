#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KNN graph calculation.

Builds a KNN graph with PCA-transformed data.
"""
import scanpy as sc
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

def main(config: dict) -> None:
    """Build a weighted KNN graph from PCA embeddings"""

    input_path = Path(config["paths"]["processed_dir"]) / "pca_result.h5ad"
    output_path = Path(config["paths"]["processed_dir"]) / "neighbor_graph.h5ad"
    n_neighbors = config["params"]["dim_reduction"]["n_neighbors"]
    n_pcs = config["params"]["dim_reduction"]["n_pcs"]

    print(f"Loading PCA-transformed data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # ---构建KNN图---
    print(f"Building KNN graph (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # ---验证图像构建---
    print(f"Distances matrix shape: {adata.obsp['distances'].shape}")
    print(f"Connectivities matrix shape: {adata.obsp['connectivities'].shape}")

    print(f"Saving graph-embedded data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    config = get_config()
    main(config)
    print("Neighbor graph construction completed successfully.")