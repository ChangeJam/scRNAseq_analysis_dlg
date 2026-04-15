import scanpy as sc
import os

TARGET_SUM  = 1e4

def normalize_data(input_path: str, output_path: str) -> None:
    print(f"Loading data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    print("Applying log1p transformation to stablize variance...")
    sc.pp.log1p(adata)

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