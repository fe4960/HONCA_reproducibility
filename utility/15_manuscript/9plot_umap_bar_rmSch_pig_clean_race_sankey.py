import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import rcParams
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/figures/"
########
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean"
bname="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_rmUnk_clean"
adata=sc.read(f"{outdir}/{bname}.h5ad")


rcParams["figure.figsize"] = (5,5)
rcParams.update({'font.size': 14})


meta=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid", sep=",", header=0)
meta=meta.set_index("sampleid")



#def generate_txt(outdir, bname):
adata=sc.read(f'{outdir}/{bname}.h5ad')
adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]
#obs=adata.obs.copy()
dt={"Asian-Indian": "Asian"}
mask=adata.obs["source"] == "Sanes"
adata.obs.loc[mask, "race"] = adata.obs.loc[mask, "sampleid"].map(meta["race"])
adata.obs["race"]=adata.obs["race"].astype("str").replace(dt)

adata.obs[["nCount_RNA","majorclass1","celltype"]].groupby(["majorclass1","celltype"]).count().to_csv(f"{outdir}/{bname}_sankey.csv",sep="\t")

sc.pl.umap(adata, color="race", frameon=False, palette="tab20", save=f'_{bname}_race_corrected.png')


