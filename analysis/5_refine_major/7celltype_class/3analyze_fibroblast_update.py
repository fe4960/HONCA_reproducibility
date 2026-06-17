import pandas as pd
import numpy as np
import scanpy as sc
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
bname0="Fibroblast_hvg_10000_fl_seurat_res_1_sd_7_clean_final_update"

bname="Fibroblast_hvg_10000_fl_seurat_res_1_sd_7_clean_final_update"



adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/{bname0}_scvi_trg.h5ad")
obs=pd.read_csv(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/{bname}.obs.gz",header=0,index_col=0)
adata.obs["leiden"]=obs.loc[adata.obs.index,"leiden"].astype("str")
df=adata.obs.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=0).values[np.newaxis,:]
conf_mat.to_csv(f'{outdir}/{bname}_conf_mat_table.csv',header=True, index=True)


#row_max=conf_mat.max(axis=1)
#dominate_name=row_max[row_max > 0.7].index
#print(dominate_name)
#conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_conf_mat_table_max.csv")
#sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"
#sc.pl.embedding(adata, basis="X_umap", color=["celltype"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_celltype_onData.png', palette="tab20")

#data=pd.read_csv(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/{bname}_conf_mat_table.csv",header=0)


df=adata.obs.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]

conf_mat.to_csv(f'{outdir}/{bname}_conf_mat_table_verticle.csv',header=True, index=True)
