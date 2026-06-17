import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

#rcParams["figure.figsize"] = (5,5)
dir="/dfs3b/ruic20_lab/junw42"

#sc.settings.figdir=f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/figures"

adata=sc.read(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge.h5ad")

obs=pd.read_csv(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg_10000_epoch_none_seurat_v3_rs_0.4_res_0.4.obs.gz",header=0,index_col=0)
adata.obs["leiden"]=obs.loc[adata.obs.index,"leiden"]

adata1=adata[adata.obs["leiden"]!=3]
adata1.write(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_rawcount_rm_res0.4_clu3.h5ad")

adata2=sc.read("HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg10000_epochnone_seurat_v3_rs_0.4_rm3_none_scvi_trg.h5ad")
adata2.obs.to_csv("HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg10000_epochnone_seurat_v3_rs_0.4_rm3_none_scvi_trg.obs.gz")
