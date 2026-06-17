import scanpy as sc
adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_kp_cluster_2_seed_7_scvi_trg.h5ad")
import pandas as pd
obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_kp_cluster_2_seed_7.obs.gz",header=0,index_col=0)
adata.obs["leiden1"]=obs.loc[adata.obs.index,"leiden"].astype("str")
obs.loc[obs["leiden"].isin([2,4]),].to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_kp_cluster_2_seed_7_kp_only_2_4.obs.gz")

clu24=obs.loc[obs["leiden"].isin([2,4]),]

adata2=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_15_16_seed12345_rs0.8_rm5_9_seed_7_scvi_trg.h5ad")

adata2.obs.loc[~adata2.obs.index.isin(clu24.index)].to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_15_16_seed12345_rs0.8_rm5_9_seed_7_rm_retonh_clu2_4.obs.gz")


clu2=obs.loc[obs["leiden"]==2,]
clu4=adata1.obs.loc[adata1.obs["leiden"]=="4",]
adata2.obs.loc[~(adata2.obs.index.isin(clu2.index) | adata2.obs.index.isin(clu4.index))].to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_15_16_seed12345_rs0.8_rm5_9_seed_7_rm_retonh_clu2_clu4.obs.gz")

