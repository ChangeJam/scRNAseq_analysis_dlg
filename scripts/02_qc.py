import scanpy as sc
import numpy as np

#Load data
adata = sc.read_h5ad("/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data/raw_data.h5ad")

#The Known mitochondrial gene IDs
mt_genes = ['FBgn0015553', 'FBgn0015558', 'FBgn0015561', 'FBgn0015550']

#Create a new list in adata.var, note which is mitochondrial genes
adata.var['mt'] = adata.var_names.isin(mt_genes)

# Calculate QC metrics. qc_vars=['mt'] calculates the percentage of mitochondrial counts, 
# which is stored in adata.obs['pct_counts_mt']
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

#Set the default directory for saving figures 
sc.settings.figdir = "/data2/ChangeJam/scRNA_seq_scrib_dlg/figures/qc_plots"

sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save='_before_filter.pdf'
)

sc.pl.scatter(
    adata, 
    x='total_counts', 
    y='n_genes_by_counts', 
    color='pct_counts_mt',
    save='_before_filter.pdf'
)

def mad_outlier(adata, metric, n_thresh=5):
    data = adata.obs[metric]
    #Calculate the median
    med = np.median(data)
    #Calculate the Median Absolute Deviation (MAD)
    mad = np.median(np.abs(data - med)) * 1.4826

    #Define the upper and lower bounds
    upper_bound = med + n_thresh * mad
    lower_bound = med - n_thresh * mad

    #Return a Boolean array marking cells that exceed the bounds as True
    return (data > upper_bound) | (data < lower_bound)

print(f"n_obs(filter before): {adata.n_obs}")

#Mark the outliers using the MAD method (here, n_thresh=5 is relatively loose)
outlier_genes = mad_outlier(adata, 'n_genes_by_counts', n_thresh=5)
outlier_counts = mad_outlier(adata, 'total_counts', n_thresh=5)
outlier_mt = mad_outlier(adata, 'pct_counts_mt', n_thresh=5)

#Use bitwise NOT (~) to keep only normal (non-outlier) cells
adata = adata[~(outlier_genes | outlier_counts | outlier_mt), :].copy()

print(f"n_obs(filtered): {adata.n_obs}")

adata.write_h5ad("/data2/ChangeJam/scRNA_seq_scrib_dlg/raw_data/filtered_data.h5ad")
print("QC done, data has been saved!")