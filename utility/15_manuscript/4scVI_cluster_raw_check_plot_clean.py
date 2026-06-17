import scanpy as sc
import pandas as pd
import numpy as np

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)
rcParams['font.family'] = 'sans-serif'
rcParams["font.sans-serif"] = ["Nimbus Sans"]

cellnum=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new.csv",header=0,index_col=0)

cellnum1=cellnum.loc[-cellnum["major class"].isin(["BC","Rod","Cone","HC","AC","RGC","Adipocyte","Pigmented_cell","RPE","Schwann_cell"]),]


bname="major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/"
sc.settings.figdir = f"{outdir}/figures"
adata=sc.read(f"{outdir}/major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet_scvi_trg.h5ad")

adata=adata[~adata.obs["majorclass"].isin(["BC","Rod","Cone","HC","AC","RGC","Adipocyte","Pigmented_cell","RPE","Schwann_cell"])]

group=cellnum1["major class"]  ##["Oligodendrocyte","Astrocyte","Fibroblast",]
sc.pl.umap(adata,color="majorclass",save=f"_{bname}_clean_majorclass.png", frameon=False, palette="tab20")


group=["Astrocyte","Endothelial_cell","Fibroblast","MG","Macrophage","Melanocyte","Microglia","Mural_cell","Oligodendrocyte","Oligodendrocyte_precursor_cell","T_cell","NK_cell","Mast_cell","Dendritic_cell","B_cell"]
mk=["GFAP","AC092957.1","MECOM","VWF","BICC1","GABRG3","HKDC1","F13A1","PAX3","KCNQ3","ADAM28","MYH11","CARMN","RNF220","CTNNA3","NXPH1","CSMD1","BCL11B","PYHIN1","CD247","KIT","SKAP1","BANK1"]
sc.pl.dotplot(adata,mk,groupby="majorclass",categories_order=group, save=f"_{bname}_clean_majorclass.png")

