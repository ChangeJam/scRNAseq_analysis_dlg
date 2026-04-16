import scanpy as sc
import os

# Standard number of highly variable genes to retain
N_TOP_GENES = 2000

def find_and_subset_hvgs(input_path: str, output_path: str) -> None:
    """
    Docstring: Identify HVGs based on normalized data and subset the AnnData object
    """

    # Load data
    print(f"Loading normalized data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Calculate HVGs
    print(f"Calculating top {N_TOP_GENES} highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        flavor='seurat',
        subset=False  # Do not slice the matrix immediately
    )

    # Subset matrix using the boolean mask on genes
    print("Subsetting AnnData to retain only HVGs...")
    adata = adata[:, adata.var.highly_variable]

    print(f"Remaining genes after subsetting: {adata.n_vars}")

    # Output the result
    print(f"Saving HVG-subset data to: {output_path}")
    adata.write_h5ad(output_path)


if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "normalized_data.h5ad"
    OUTPUT_FILE = "hvg_subset.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)
    output_path = os.path.join(DATA_DIR, OUTPUT_FILE)

    find_and_subset_hvgs(input_path, output_path)
    print("HVG selection and subsetting completed successfully.")