import scanpy as sc
import os


N_COMPS = 50

def run_pca(input_path: str, output_path: str) -> None:
    """
    Execute PCA and save the resulting AnnData object.
    """

    print(f"Loading HVG-subset data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    print(f"Running PCA to compute {N_COMPS} principal components...")
    sc.pp.pca(adata, n_comps=N_COMPS)

    print(f"PCA output matrix shape (Cells x PCs): {adata.obsm['X_pca'].shape}")
    print(f"PC loadings matrix shape (Genes x PCs): {adata.varm['PCs'].shape}")

    print(f"Saving PCA-transformed data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "hvg_subset.h5ad"
    OUTPUT_FILE = "pca_result.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)
    output_path = os.path.join(DATA_DIR, OUTPUT_FILE)

    run_pca(input_path, output_path)
    print("PCA computation completed successfully.")