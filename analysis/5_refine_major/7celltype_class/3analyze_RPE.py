import pandas as pd
import numpy as np
import scanpy as sc
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/"
#bname="RPE_seed_7_res_0.4_clean_other_celltype"
bname="RPE_seed_7_res_0.4_clean_other_celltype_rm"

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/major_seed_7_res_0.3_clean_other_scvi_trg.h5ad")


adata1=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/RPE_hvg2000_epochnone_seurat_v3_rs_0.4_clean_rm_scvi_trg.h5ad")

idx0=adata[adata.obs["majorclass"]=="RPE"].obs.index.difference(adata1.obs.index)

adata=adata[~adata.obs.index.isin(idx0)]

idx=adata.obs.index.intersection(adata1.obs.index)
adata.obs["celltype"]=adata.obs["celltype"].cat.add_categories("4")
adata.obs.loc[idx,"celltype"]=adata1.obs.loc[idx,"leiden"].astype(str)


adata.obs.loc[idx,"celltype"]=adata1.obs.loc[idx,"leiden"]

df=adata.obs.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=0).values[np.newaxis,:]
conf_mat.to_csv(f'{outdir}/{bname}_conf_mat_table.csv',header=True, index=True)


row_max=conf_mat.max(axis=1)
dominate_name=row_max[row_max > 0.7].index
print(dominate_name)
conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_conf_mat_table_max.csv")
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/figures"
sc.pl.embedding(adata, basis="X_umap", color=["celltype"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_celltype_onData.png', palette="tab20")

data=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/RPE_seed_7_res_0.4_clean_other_conf_mat_table.csv",header=0)
df=data.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
#outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/"
#conf_mat.to_csv(f'{outdir}/{bname}_celltype_conf_mat_table.csv',header=True, index=True)
