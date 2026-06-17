import scanpy as sc
import pandas as pd
adata=sc.read("HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg2000_epochnone_seurat_v3_rs_0.4_clean_rm4_sb_scvi_trg.h5ad")
res=[0.7,0.8,0.9, 0.1]
for rs in res:
    obs=pd.read_csv("")
