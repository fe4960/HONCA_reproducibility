import scanpy as sc
import pandas as pd
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad")
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures/"

mk=["PIEZO2","PIEZO1", "BMP4", "BMP5", "BMP7"]

sc.pl.umap(adata, color=mk, ncols=1,  save= "Astrocyte_subclass_new5_clean_PIEZO_BMP.png" )


