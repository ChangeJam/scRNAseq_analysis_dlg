#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute UMAP embeddings from the KNN graph.
Projects 50-dimensional PCA space into 2 dimensions for visualization.
"""

import scanpy as sc
import os

def run_umap(input_path: str, output_path: str) -> None:
    """Execute UMAP dimensionality reduction and save the result."""

    print(f"Loading clustered data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Run UMAP
    print("Running UMAP...")
    sc.tl.umap(adata)

    # Verify output
    print(f"UMAP embedding shape: {adata.obsm['X_umap'].shape}")

    print(f"Saving UMAP-transformed data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "clustered.h5ad"
    OUTPUT_FILE = "umap_result.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)
    output_path = os.path.join(DATA_DIR, OUTPUT_FILE)

    run_umap(input_path, output_path)
    print("UMAP computation completed successfully.")