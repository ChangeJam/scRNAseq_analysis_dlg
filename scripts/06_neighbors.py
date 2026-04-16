import scanpy as sc
import os

# Each cell connects to its 15 nearest neighbors
N_NEIGHBORS = 15

# Use the top 30 PCs for distance calculation
N_PCS = 30

def build_neighbor_graph(input_path: str, output_path: str) -> None:
    """Build a weighted KNN graph from PCA embeddings"""

    print(f"Loading PCA-transformed data from: {input_path}")
    adata = sc.read_h5ad(input_path)

    # Build the neighbor graph
    print(f"Building KNN graph (n_neighbors={N_NEIGHBORS}, n_pcs={N_PCS})...")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)

    # Verify graph construction
    print(f"Distances matrix shape: {adata.obsp['distances'].shape}")
    print(f"Connectivities matrix shape: {adata.obsp['connectivities'].shape}")

    print(f"Saving graph-embedded data to: {output_path}")
    adata.write_h5ad(output_path)

if __name__ == "__main__":
    DATA_DIR = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data"
    INPUT_FILE = "pca_result.h5ad"
    OUTPUT_FILE = "neighbor_graph.h5ad"

    input_path = os.path.join(DATA_DIR, INPUT_FILE)
    output_path = os.path.join(DATA_DIR, OUTPUT_FILE)

    build_neighbor_graph(input_path, output_path)
    print("Neighbor graph construction completed successfully.")