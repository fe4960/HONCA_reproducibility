import scanpy as sc
import pandas as pd
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"
mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/LC",header=None)
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass.h5ad")
sc.pl.dotplot(adata,mk[0],groupby="subclass", save="Fibroblast_subclass_LC.png")
