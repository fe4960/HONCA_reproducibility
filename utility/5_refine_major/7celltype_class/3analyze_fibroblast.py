import pandas as pd
import numpy as np
import scanpy as sc
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
#bname="Fibroblast_seed_7_res_0.4_clean_other_celltype"
#bname="Fibroblast_seed_7_res_0.4_clean_other_celltype_rm"

bname="Fibroblast_hvg_10000_fl_seurat_v3_res_0.4_sd_7_clean_final"


#bname="Fibroblast_hvg_2000_fl_seurat_v3_res_0.4_sd_7_clean_other"
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_seed_7_res_0.4_clean_other_scvi_trg.h5ad")
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_seed_7_res_0.4_clean_other_clean_scvi_trg.h5ad")
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg_2000_fl_seurat_v3_res_0.4_sd_7_clean_other_scvi_trg.h5ad")

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg_10000_fl_seurat_v3_res_0.4_sd_7_clean_final_scvi_trg.h5ad")

df=adata.obs.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=0).values[np.newaxis,:]
conf_mat.to_csv(f'{outdir}/{bname}_conf_mat_table.csv',header=True, index=True)


row_max=conf_mat.max(axis=1)
dominate_name=row_max[row_max > 0.7].index
print(dominate_name)
conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_conf_mat_table_max.csv")
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"
sc.pl.embedding(adata, basis="X_umap", color=["celltype"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_celltype_onData.png', palette="tab20")

data=pd.read_csv(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/{bname}_conf_mat_table.csv",header=0)


#data=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_seed_7_res_0.4_clean_other_conf_mat_table.csv",header=0)
df=data.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
#outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
#conf_mat.to_csv(f'{outdir}/{bname}_celltype_conf_mat_table.csv',header=True, index=True)
