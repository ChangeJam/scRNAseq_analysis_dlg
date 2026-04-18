#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load raw 10x Genomics data and convert to h5ad format.

Reads the matrix.mtx.gz / barcode.tsv.gz / features.tsv.gz triplet
produced by CellRanger and saves a single AnnData (.h5ad) file.
"""

import scanpy as sc
from pathlib import Path
import sys

# --- 使 Python 能找到 src/ 包 ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.config import get_config


def load_data(config: dict) -> sc.AnnData:
    """Read 10x MTX data and save as h5ad"""

    raw_dir = Path(config["paths"]["raw_data_dir"])
    processed_dir = Path(config["paths"]["processed_dir"])

    # 确保输出目录存在
    processed_dir.mkdir(parents=True, exist_ok=True)

    output_path = processed_dir / "raw_data.h5ad"

    # 读取10x三个文件
    print(f"Reading 10x data from: {raw_dir}")
    adata = sc.read_10x_mtx(raw_dir)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # 保存为h5ad
    print(f"Saving to: {output_path}")
    adata.write_h5ad(output_path)

    return adata

if __name__ == "__main__":
    config = get_config()
    load_data(config)
    print("Step 01 completed: raw data loaded.")