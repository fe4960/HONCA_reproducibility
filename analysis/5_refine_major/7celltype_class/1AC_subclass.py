import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/AC/clean/figures"

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/AC/"
mk="sanes_mk"
adata_query=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/query_query_latent.h5ad")

adata_query1=sc.read_h5ad(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/AC/clean/AC_res_1_clean_res_1.0.h5ad")


celltype="AC"
idx=adata_query[adata_query.obs["scANVI_prediction_max_probability"]>0.9].obs.index.intersection(adata_query1.obs.index)

adata_query1=adata_query1[idx]
#adata_query1.obs.loc[idx,"subclass"]=adata_query.obs.loc[idx,"scANVI_predictions"]
adata_query1=adata_query1[-adata_query1.obs["leiden"].isin(["43","42"])] #43 is contaminated by OLIG, 42 is wired.
adata_query1.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/AC/clean/AC_res_1_clean_rm_res_1.0.h5ad")
#######

