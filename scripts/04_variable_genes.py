#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Highly Variable Genes (HVG) selection 

Identifies HVGs based on normalized data, subsets the 
AnnData object to retain only these genes, and save the result.
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
    Identify HVGs based on normalized data and subset the AnnData object.
    """

    input_path = Path(config["paths"]["processed_dir"]) / "normalized_data.h5ad"
    output_path = Path(config["paths"]["processed_dir"]) / "hvg_subset.h5ad"
    n_top_genes = config["params"]["hvg"]["n_top_genes"]
    flavor = config["params"]["hvg"]["flavor"]

    # Load data
    print(f"Loading normalized data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Calculate HVGs
    print(f"Calculating top {n_top_genes} highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        subset=False  # Do not slice the matrix immediately
    )

    # Subset matrix using the boolean mask on genes
    print("Subsetting AnnData to retain only HVGs...")
    adata = adata[:, adata.var.highly_variable]

    print(f"Remaining genes after subsetting: {adata.n_vars}")

    # Save subsetted data
    print(f"Saving HVG-subset data to: {output_path}")
    adata.write_h5ad(output_path)


if __name__ == "__main__":
    config = get_config()
    main(config)
    print("Highly variable genes completed successfully!")