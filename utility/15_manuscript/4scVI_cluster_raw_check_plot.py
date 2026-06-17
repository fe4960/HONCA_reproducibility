import scanpy as sc
import pandas as pd
import numpy as np

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

bname="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE"
#bname="major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/"
sc.settings.figdir = f"{outdir}/figures"
#adata_query=sc.read(f"{outdir}/major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet_scvi_trg.h5ad")

adata_query=sc.read(f"{outdir}/major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_scvi_trg.h5ad")


adata_query=adata_query[-adata_query.obs["majorclass1"].isin(["Adipocyte", "Schwann_cell", "Pigmented_cell"])]
sc.pl.embedding(adata_query, basis="X_umap", color=["majorclass"], ncols=1,frameon=False,save=f'_{bname}_majorclass_onData.png', palette="tab20")

dt={"T_cell": "Immune_cell", "NK_cell": "Immune_cell","Mast_cell":"Immune_cell", "Dendritic_cell": "Immune_cell", "B_cell": "Immune_cell"}

adata_query.obs["majorclass_wRet"]=adata_query.obs["majorclass"].replace(dt)

sc.pl.embedding(adata_query, basis="X_umap", color=["majorclass_wRet"], ncols=1,frameon=False,save=f'_{bname}_majorclass_wRet_onData.png', palette="tab20")

mk=["AQP7","PLIN1"]

sc.pl.dotplot(adata_query, mk, groupby="leiden",save=f'_{bname}_majorclass_leiden_adipocyte.png')

mk=["CCBE1","CPNE5"]

sc.pl.dotplot(adata_query, mk, groupby="leiden",save=f'_{bname}_majorclass_leiden_pigmented_cell.png')


ct=pd.read_csv(f"/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/9_cellxgene/5leiden_celltype",sep="\t",header=None)
dt=ct.set_index(0).to_dict()[1]
#bname="ONONH_hvg2k_epoch20_new"
adata_query.obs["majorclass_woRet"]=adata_query.obs["leiden"].astype(int).replace(dt,regex=True)

#adata_query1=sc.read_h5ad(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/cellclass/MG/MG_hvg2k_epoch20_scvi.h5ad")
adata_query1=sc.read_h5ad(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/major_hvg10000_epochnone_seurat_v3_rs_1_sb_kp_cl_16_MG_scvi_trg.h5ad")

idx=adata_query1[adata_query1.obs["leiden"].isin(["4"])].obs.index.intersection(adata_query.obs.index)


#idx=adata_query1[adata_query1.obs["leiden"].isin(["4","9"])].obs.index.intersection(adata_query.obs.index)
adata_query.obs.loc[idx,"majorclass_woRet"]="Astrocyte"


sc.pl.embedding(adata_query, basis="X_umap", color=["majorclass_woRet"], ncols=1,frameon=False,save=f'_{bname}_majorclass_woRet_onData.png', palette="tab20")

df = adata_query.obs.groupby(["majorclass_woRet", "majorclass_wRet"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)
norm_df.sort_index(inplace=True)
norm_df.sort_index(inplace=True,axis=1)
norm_df = norm_df.reindex(sorted(norm_df.columns), axis=1)

norm_df.to_csv(f'{outdir}/{bname}_wRet_woRet_conf_mat_table.csv',header=True, index=True)

import matplotlib.pyplot as plt

fig=plt.figure(figsize=(15, 15))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("majorclass_wRet")
plt.ylabel("majorclass_woRet")
fig.savefig(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/figures/{bname}_wRet_woRet_conf.png')



