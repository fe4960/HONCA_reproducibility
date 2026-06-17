import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
celltype="Microglia"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
mk="sanes_mk"
adata=sc.read_h5ad(f"{outdir}/{celltype}_subclass_sb_clean.h5ad")

bname=f"{celltype}_subclass_sb_clean"



sc.settings.figdir = f"{outdir}/figures"


sc.pl.embedding(adata, basis="X_umap", color=["subclass"],  ncols=1,frameon=False,save=f'{bname}.png', palette="tab20")


