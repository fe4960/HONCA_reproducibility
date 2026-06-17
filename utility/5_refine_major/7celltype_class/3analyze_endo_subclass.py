adata1=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_subclass_sb.h5ad")
import pandas as pd
mk="NF2"
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/figures"
sc.pl.dotplot(adata1,mk,groupby="subclass", save="Endothelial_cell_subclass_sb_meningothelial.png")
