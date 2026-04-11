import scanpy as sc

# 1. Define Path
mtx_dir = "/data2/ChangeJam/scRNA_seq_scrib_CellCompetition/1.Expression/expressions/dlg"
output_h5ad = "/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data/raw_data.h5ad"

# 2. Read mtx format data
print("开始读取数据...")
adata = sc.read_10x_mtx(mtx_dir)
print(f"读取成功！共 {adata.n_obs} 个细胞，{adata.n_vars} 个基因")

# 3. Save in h5ad format
print(f"正在保存到 {output_h5ad} ...")
adata.write_h5ad(output_h5ad)
print("保存成功！")