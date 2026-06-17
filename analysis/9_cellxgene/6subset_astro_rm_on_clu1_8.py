import scanpy as sc
import pandas as pd
celltype="Astrocyte"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
adata=sc.read_h5ad(f"{outdir}/{celltype}_subclass_new2.h5ad")
obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2_on_hvg10000_epochnone_seurat_v3_rs_0.4_clean_sb.obs.gz",header=0,index_col=0)
#adata=adata[adata.obs["subclass"].isin(["Astro_retina_NLGN1+","Astro_retina","Astro_ONH"])]
#adata.write(f"{outdir}/{celltype}_subclass_new2_ret_onh.h5ad" )
idx=obs.loc[obs["leiden"].isin([1,8])].index

adata=adata[~adata.obs.index.isin(idx)]

adata.write(f"{outdir}/{celltype}_subclass_new2_rm_clu1_8.h5ad")
