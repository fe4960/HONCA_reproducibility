import scanpy as sc

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg.h5ad")

adata=adata[-adata.obs["leiden"].isin(["9","10"])]

import pandas as pd
obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.2_clean_sb_rm0_12_15_16_merge_seed_7_kp_0_res_0.5_scvi_trg.h5ad.obs.gz",header=0,index_col=0)
obs["leiden"]="on_"+obs["leiden"].astype("str")

adata1=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_none_rm2_seed_7_scvi_trg.h5ad")

adata1=adata1[-adata1.obs["leiden"].isin(["4"])]
obs1=adata1.obs.copy()
obs1["leiden"]="ret_"+obs1["leiden"].astype("str")

idx=obs1.index.difference(obs.index)

obs_full=pd.concat([obs,obs1.loc[idx,]])


adata.obs["leiden1"]=adata.obs["leiden"].astype("str")

idx=adata.obs.index.intersection(obs_full.index)

adata.obs.loc[idx,"leiden1"]=obs_full.loc[idx,"leiden"].astype("str")


adata=adata[-adata.obs["leiden1"].isin(["1","2","3","4","5","6"])].copy()

#index=adata.obs.index.intersection(obs_full.index)



#index=index.unique()

#adata.obs.loc[index,"leiden1"]=obs_full.loc[index,"leiden"].astype("str")

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"

celltype="Astrocyte"

adata_full=sc.read_h5ad(f"{outdir}/{celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
#bname=f"{celltype}_subclass_new5_clean"

#adata_full.obs["subclass"]="Astro"
idx=adata.obs.index
#adata=adata[index].copy()
adata_full=adata_full[idx].copy()
#sc.settings.figdir = f"{outdir}/figures"


del adata_full.obsm["X_scVI"]
adata_full.obsm["X_scVI"]=adata.obsm["X_scVI"]
del adata_full.obsm["X_umap"]
adata_full.obsm["X_umap"]=adata.obsm["X_umap"]
del adata_full.uns["umap"]
adata_full.uns["umap"]=adata.uns["umap"]
del adata_full.obs["leiden"]


adata_full.obs["leiden1"]=adata.obs["leiden1"].astype("str")
bname="Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_full_clean_update"
adata_full.write(f"{outdir}/{bname}.h5ad")
