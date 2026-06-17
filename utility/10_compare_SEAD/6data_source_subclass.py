import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

rcParams["figure.figsize"] = (5,5)
dir="/dfs3b/ruic20_lab/junw42"

sc.settings.figdir=f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/figures"

adata=sc.read(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg_10000_epoch_none_seurat_v3_rs_1_scvi_trg.h5ad")
sc.pl.umap(adata,color="subclass",save="_Oligo_GSE267301_merge_subclass.png", palette=dt,size=1, alpha=0.5)
