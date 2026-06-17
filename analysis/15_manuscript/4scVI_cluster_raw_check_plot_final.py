import scanpy as sc
import pandas as pd
import numpy as np

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

bname="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE"
#bname="major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/"
sc.settings.figdir = f"{outdir}/figures"
#adata_query=sc.read(f"{outdir}/major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet_scvi_trg.h5ad")

adata_query=sc.read(f"{outdir}/major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_scvi_trg.h5ad")


adata_query=adata_query[(-adata_query.obs["majorclass1"].isin(["Adipocyte", "Schwann_cell", "Pigmented_cell"])) & (adata_query.obs["celltype"]!="Unknown?")]
sc.pl.embedding(adata_query, basis="X_umap", color=["majorclass1"],   ncols=1,frameon=False,save=f'_{bname}_majorclass_onData.png', palette="tab20")


sc.pl.embedding(adata_query, basis="X_umap", color=["celltype"], ncols=1,frameon=False,save=f'_{bname}_celltype_onData.png', palette="tab20")
