#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Perform Leiden clustering on the KNN graph.
Assigns each cell to a cluster label based on graph connectivity density.
"""

import scanpy as sc
import os

# Controls clustering granularity by resolution = 1.0
RESOLUTION = 1.0

def run_clustering(input_path: str, output_path: str) -> None:
    """Execute Leiden clustering and save the result."""

    print(f"Loading neighbor-graph data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Run Leiden clustering
    print(f"Running Leiden clustering (resolution={RESOLUTION})...")
    sc.tl.leiden(adata, resolution=RESOLUTION)

    # Show cluster summary
    print("Cluster distribution:")
    print(adata.obs['leiden'].value_counts().sort_index())

    print(f"Saving clustered data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "neighbor_graph.h5ad"
    OUTPUT_FILE = "clustered.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)
    output_path = os.path.join(DATA_DIR, OUTPUT_FILE)

    run_clustering(input_path, output_path)
    print("Leiden clustering completed successfully.")