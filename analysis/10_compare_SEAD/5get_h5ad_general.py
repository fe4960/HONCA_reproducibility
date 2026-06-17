import scanpy as sc
import anndata as ad
import sys
import pandas as pd

#donor={'SRR29007591': 'ALSP',
# 'SRR29007592': 'ALSP',
# 'SRR29007593': 'ALSP',
# 'SRR29007594': 'ALSP',
# 'SRR29007595': 'ALSP',
# 'SRR30562939': 'AD',
# 'SRR30562940': 'Ctr',
# 'SRR30562941': 'Ctr',
# 'SRR30562942': 'Ctr',
# 'SRR30562943': 'Ctr',
# 'SRR30562944': 'Ctr',
# 'SRR30562945': 'AD'}


adata0=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/CSF1R/GSE267301/GSE267301_hvg_10000_epoch_none_seurat_v3_rs_1_scvi_trg.h5ad")
adata_full0=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/CSF1R/GSE267301.h5ad")

adata0=adata0[adata0.obs["sampleid"].isin(['SRR30562940','SRR30562941','SRR30562942','SRR30562943','SRR30562944'])]
adata_full0=adata_full0[adata_full0.obs["sampleid"].isin(['SRR30562940','SRR30562941','SRR30562942','SRR30562943','SRR30562944'])]


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/GSE167494/GSE167494_hvg_10000_epoch_none_seurat_v3_rs_1_sb_scvi_trg.h5ad")
adata_full=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/GSE167494/GSE167494.h5ad")

cell=sys.argv[1]
fn0=sys.argv[2]
fn1=sys.argv[3]
subclass=sys.argv[4]

list0=pd.read_csv(fn0,header=None)
leiden0=list0[0].values.astype("str")

list1=pd.read_csv(fn1,header=None)
leiden1=list1[0].values.astype("str")

astro0=adata0[adata0.obs["leiden"].isin(leiden0)].obs.index

adata_astro0=adata_full0[astro0].copy()
adata_astro0.obs["majorclass"]=cell
adata_astro0.obs["subclass"]=f"{cell}_GSE267301"



astro=adata[adata.obs["leiden"].isin(leiden1)].obs.index

adata_astro=adata_full[astro].copy()
adata_astro.obs["majorclass"]=cell
adata_astro.obs["subclass"]=f"{cell}_GSE167494"


h5ad=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{cell}/clean/{subclass}.h5ad"

adata1=sc.read(h5ad)
adata1.X=adata1.layers["counts"]

adata_full1=ad.concat([adata_astro,adata1])


adata_full1.write(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/{cell}_PFC_ONONH.h5ad")


