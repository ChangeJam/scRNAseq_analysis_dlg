import scanpy as sc
import os

# Standard scaling factor for single-cell library size normalization
TARGET_SUM = 1e4

def normalize_data(input_path: str, output_path: str) -> None:
    # Load filtered data, perform normalization and log transformation, then save

    # Load data
    print(f"Loading data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Normalize total counts per cell
    # Scale raw counts to a common target sum across all cells
    print(f"Normalizing total counts to {TARGET_SUM} per cell...")
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)

    # Apply log(1+x) to handle zero inflation and stabilize variance
    print("Applying log1p transformation to stabilize variance...")
    sc.pp.log1p(adata)

    # Save results
    print(f"Saving normalized data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "filtered_data.h5ad"
    OUTPUT_FILE = "normalized_data.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)
    output_path = os.path.join(DATA_DIR, OUTPUT_FILE)

    normalize_data(input_path, output_path)
    print("Normalization completed successfully.")