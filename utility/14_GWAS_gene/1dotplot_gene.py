import scanpy as sc

import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

rcParams["figure.figsize"] = (5,5)


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_normcount.h5ad")

dir="/dfs3b/ruic20_lab/junw42"
bname=f"major_ON_ONH_clean_new"
sc.settings.figdir = f"{dir}/HCA_ON/data/5_refine_major/scvi/major/GWAS/figures"

sc.pl.violin(adata,["CLIC1"],groupby="majorclass",save=f"_{bname}_CLIC1.png", rotation=90)
sc.pl.dotplot(adata,["CLIC1"],groupby="majorclass",save=f"_{bname}_CLIC1.png")


