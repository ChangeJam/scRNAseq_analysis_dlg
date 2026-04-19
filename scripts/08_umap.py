#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute UMAP embeddings from the KNN graph.

Projects PCA embeddings into 2 dimensions for visualization.
"""

import scanpy as sc
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

def main(config: dict) -> None:
    """Execute UMAP dimensionality reduction and save the result."""

    input_path = Path(config["paths"]["processed_dir"]) / "clustered.h5ad"
    output_path = Path(config["paths"]["processed_dir"]) / "umap_result.h5ad"

    print(f"Loading clustered data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # ---运行UMAP---
    print("Running UMAP...")
    sc.tl.umap(adata)

    # ---验证输出---
    print(f"UMAP embedding shape: {adata.obsm['X_umap'].shape}")

    print(f"Saving UMAP-transformed data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    config = get_config()
    main(config)
    print("UMAP computation completed successfully.")