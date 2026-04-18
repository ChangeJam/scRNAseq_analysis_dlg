#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Normalization: Perform normalization and log transformation

Reads filtered data, runs normalization and saves transformed data.
"""

import scanpy as sc
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config

# Standard scaling factor for single-cell library size normalization
# TARGET_SUM = 1e4

def main(config: dict) -> None:
    """Load filtered data, perform normalization and log transformation, then save"""

    input_path = Path(config["paths"]["processed_dir"]) / "filtered_data.h5ad"
    output_path = Path(config["paths"]["processed_dir"]) / "normalized_data.h5ad"
    target_sum = config["params"]["normalization"]["target_sum"]

    # ---导入数据---
    print(f"Loading data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # ---标准化每一个细胞---
    print(f"Normalizing total counts to {target_sum} per cell")
    sc.pp.normalize_total(adata, target_sum = target_sum)

    # ---log转换---
    print("Applying log1p transformation to stabilize variance...")
    sc.pp.log1p(adata)

    # ---结果存储---
    print(f"Saving normalized data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    config = get_config()
    main(config)
    print("Normalization completed successfully.")