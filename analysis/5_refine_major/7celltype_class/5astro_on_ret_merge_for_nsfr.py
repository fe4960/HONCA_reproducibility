import scanpy as sc
from matplotlib import rcParams
import pandas as pd
import matplotlib.pyplot as plt
dir1="/dfs3b/ruic20_lab/junw42"
celltype="Astrocyte"
#obs=pd.read_csv(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_none_rm4_9_rm2_seed_7.obs.gz",header=0,index_col=0)
adata=sc.read(f"{dir1}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/{celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
adata1=sc.read(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg_update.h5ad")
adata=adata[adata1.obs.index]

del adata.obsm["X_scVI"]
adata.obsm["X_scVI"]=adata1.obsm["X_scVI"]
del adata.obsm["X_umap"]
adata.obsm["X_umap"]=adata1.obsm["X_umap"]
del adata.uns["umap"]
adata.uns["umap"]=adata1.uns["umap"]
del adata.obs["leiden"]
adata.obs["leiden"]=adata1.obs["leiden1"]

adata.write(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg_update_full.h5ad")
