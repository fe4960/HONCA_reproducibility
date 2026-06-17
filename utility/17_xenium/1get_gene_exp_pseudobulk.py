import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_rawcount.h5ad")
adata=adata[-adata.obs["majorclass"].isin(["Rod","Cone","BC","AC","HC","RGC","MG","RPE","Melanocyte","Adipocyte","Schwann_cell","Pigmented_cell"])]
import numpy as np
#adata=adata[adata.obs["tissue"]=="ON"]
exp=np.ravel(adata.X.sum(axis=0, dtype=np.float64))
import pandas as pd
data=pd.DataFrame(exp)
data.index=adata.var.index
data.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/ONONH_all_rawcount_rmRet_rawgene_count.csv",sep=",")

#data.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/ONONH_all_rawcount_rmRet_rawgene_count_ON.csv",sep=",")

