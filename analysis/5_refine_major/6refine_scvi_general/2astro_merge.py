import scanpy as sc
import pandas as pd
dir1="/dfs3b/ruic20_lab/junw42"
adata=sc.read(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_none_rm4_9_seed_7_scvi_trg.h5ad")
obs2=pd.read_csv(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_none_rm4_9_seed_7.obs.gz",header=0,index_col=0)
adata.obs["leiden"]=obs2.loc[adata.obs.index,"leiden"].astype("str")
adata=adata[-adata.obs["leiden"].isin(["7","8"])]
obs=adata.obs.copy()
obs["leiden"]="on_"+obs["leiden"].astype("str")
adata=sc.read(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_none_rm2_seed_7_scvi_trg.h5ad")

adata=adata[-adata.obs["leiden"].isin(["4"])]
obs1=adata.obs.copy()
obs1["leiden"]="ret_"+obs1["leiden"].astype("str")

obs_full=pd.concat([obs,obs1])

obs_full.to_csv(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_none_rm4_9_rm2_seed_7.obs.gz")
