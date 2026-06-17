import scanpy as sc
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean"

sc.settings.figdir=f"{dir1}/figures"

#adata=sc.read(f"{dir1}/Astrocyte_subclass_new2_rm_clu1_8_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb_scvi_trg.h5ad")
adata=sc.read(f"{dir1}/Astrocyte_subclass_new2_rm_clu1_8_hvg10000_epochnone_seurat_v3_rs_0.4_seed_7_clean_sb_scvi_trg.h5ad")

import pandas as pd
obs_on=pd.read_csv(f"{dir1}/Astrocyte_subclass_new2_on_hvg10000_epochnone_seurat_v3_rs_0.4_rm1_8_clean_sb.obs.gz",header=0,index_col=0)

####adata_ret=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2_ret_onh_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb_scvi_trg.h5ad")
#####adata_ret.obs.to_csv(f"{dir1}/Astrocyte_subclass_new2_ret_onh_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb.obs.gz")

obs_ret=pd.read_csv(f"{dir1}/Astrocyte_subclass_new2_ret_onh_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb.obs.gz",header=0,index_col=0)

obs=pd.concat([obs_on,obs_ret])
adata=adata[obs.index]
adata.obs["leiden1"]="NA"
adata.obs.loc[obs_on.index,"leiden1"]="on_"+obs_on["leiden"].astype("str")
adata.obs.loc[obs_ret.index,"leiden1"]="ret_"+obs_ret["leiden"].astype("str")
sc.pl.umap(adata,color="leiden1",save="Astrocyte_subclass_new2_rm_clu1_8_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb_ret_on_leiden.png")

adata.obs["leiden2"]="comb_"+adata.obs["leiden"].astype("str")
adata.obs.loc[obs_ret.index,"leiden2"]="ret_"+obs_ret["leiden"].astype("str")
sc.pl.umap(adata,color="leiden2",save="Astrocyte_subclass_new2_rm_clu1_8_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb_ret_leiden.png")
adata.obs.to_csv(f"{dir1}/Astrocyte_subclass_new2_comb_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb.obs.gz")


obs=pd.read_csv(f"{dir1}/Astrocyte_subclass_new2_comb_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb.obs.gz",header=0,index_col=0)

adata_veloc=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_hvg_2000_raw_velocity.h5ad")

adata_veloc=adata_veloc[obs.index]
adata_veloc.obs["leiden0"]="NA"
adata_veloc.obs["leiden1"]="NA"
adata_veloc.obs["leiden2"]="NA"

adata_veloc.obs.loc[obs.index,"leiden0"]=obs["leiden"].astype(str)
adata_veloc.obs.loc[obs.index,"leiden1"]=obs["leiden1"].astype(str)
adata_veloc.obs.loc[obs.index,"leiden2"]=obs["leiden2"].astype(str)


sc.pl.violin(adata_veloc,'latent_time',groupby="leiden0",save="Astrocyte_subclass_new2_on_hvg10000_epochnone_seurat_v3_rs_0.4_rm1_8_clean_sb_leiden0_latent_time.png", rotation=90)

sc.pl.violin(adata_veloc,'latent_time',groupby="leiden1",save="Astrocyte_subclass_new2_on_hvg10000_epochnone_seurat_v3_rs_0.4_rm1_8_clean_sb_leiden1_latent_time.png", rotation=90)

sc.pl.violin(adata_veloc,'latent_time',groupby="leiden2",save="Astrocyte_subclass_new2_on_hvg10000_epochnone_seurat_v3_rs_0.4_rm1_8_clean_sb_leiden2_latent_time.png", rotation=90)
