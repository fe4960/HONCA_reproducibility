import pandas as pd
#dir1="/dfs3b/ruic20_lab/junw42/gRNA_dropoff/data/sujung/scRNA/250510_X/cellranger/crispr_analysis"
#singlet=pd.read_csv(f"{dir1}/singlet_list",header=0)

import scanpy as sc

dir="/dfs3b/ruic20_lab/junw42"

adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2_ret_onh_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb_scvi_trg.h5ad")

#mk={"other":["GAP43","LAMA2","EBF1","FOXP1"],"OPC":["PDGFRA","CSPG4","OLIG2"],"COP":["GPR17","VCAN","NEU4","BMP4","FYN","SOX6","PTPRZ1","PDCD4"],"NFOL":["TCF7L2","CEMIP2"],"MFOL":["MBP","OPALIN","CTPS1"],"MOL":["MBP","PLP1","MOG","GRM3","KLK6","PTGDS","ANXA5","HOPX"]}

#sc.pl.dotplot(adata,mk,groupby="leiden")

adata1=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_hvg_2000_raw_velocity.h5ad")

adata.obs["latent_time"]=adata1.obs.loc[adata.obs.index,"latent_time"]

adata.obs["velocity_pseudotime"]=adata1.obs.loc[adata.obs.index,"velocity_pseudotime"]

sc.pl.violin(adata,"latent_time",groupby="leiden")

adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2_rm_clu1_8_hvg10000_epochnone_seurat_v3_rs_0.4_seed_7_clean_sb_scvi_trg.h5ad")

#adata1=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_hvg_2000_raw_velocity.h5ad")

adata.obs["latent_time"]=adata1.obs.loc[adata.obs.index,"latent_time"]

adata.obs["velocity_pseudotime"]=adata1.obs.loc[adata.obs.index,"velocity_pseudotime"]

sc.pl.violin(adata,"latent_time",groupby="leiden")

#########
adata=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new3.h5ad")

obs=pd.read_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2_rm_clu1_8_hvg10000_epochnone_seurat_v3_rs_0.8_clean_sb.obs.gz",header=0,index_col=0)

adata.obs["leiden1"]=obs.loc[adata.obs.index,"leiden"].astype("str").astype("category")

mk={"retina":["GFAP", "S100B", "SOX9", "PAX2", "PAX8"],"ONH1":["CARTPT","TRPV2","TRPV1","TRPM7","PIEZO1","PIEZO2"], "ONH":["NES", "VIM", "LGALS3", "MEGF10", "ABCA1"],"ONP":["AQP4"]}

sc.pl.dotplot(adata,mk,groupby="leiden1")

######
obs=pd.read_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/twoLevel_pan_eye_Astro_seed7/Astrocyte_subclass_new3_0.4_0.2_obs.txt.gz",header=0, sep="\t",index_col=-1)
adata=adata[obs.index]
adata.obs["leiden2"]=obs.loc[adata.obs.index,"leiden_1"].astype("str").astype("category")
sc.pl.dotplot(adata,mk,groupby="leiden2")


obs=pd.read_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_hvg_2000_raw_velocity.obs.gz",sep="\t",header=0,index_col=0)

adata.obs["latent_time"]=obs.loc[adata.obs.index,"latent_time"]
adata.obs["velocity_pseudotime"]=obs.loc[adata.obs.index,"velocity_pseudotime"]
sc.pl.violin(adata,"latent_time",groupby="subclass",rotation=90)
