#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Perform Leiden clustering on the KNN graph.

Assigns each cell to a cluster label based on graph connectivity density.
"""

import scanpy as sc
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

def main(config: dict) -> None:
    """Execute Leiden clustering and save the result."""

    input_path = Path(config["paths"]["processed_dir"]) / "neighbor_graph.h5ad"
    output_path = Path(config["paths"]["processed_dir"]) / "clustered.h5ad"
    resolution = config["params"]["clustering"]["resolution"]

    print(f"Loading neighbor-graph data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # ---运行Leiden聚类---
    print(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)

    # ---聚类总结---
    print("Cluster distribution:")
    print(adata.obs['leiden'].value_counts().sort_index())

    print(f"Saving clustered data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    config = get_config()
    main(config)
    print("Leiden clustering completed successfully.")