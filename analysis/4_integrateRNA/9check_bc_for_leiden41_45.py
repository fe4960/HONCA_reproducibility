import pandas as pd
import scanpy as sc
import numpy as np
adata_query=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_hvg2k_epoch20_new_scvi.h5ad")

pnas=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/SCP_meta.csv",header=0)


samplelist=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONH_ON_public_meta.csv'
meta=pd.read_csv(samplelist,header=0)
dt=meta[["sample","sampleid"]].set_index("sample").to_dict()["sampleid"]
pnas["bc"]=pnas["NAME"].replace(dt, regex=True)
pnas=pnas.set_index("bc")
idx=pnas.index.intersection(adata_query.obs.index)

len(pnas.loc[pnas["biosample_id"]=="Hu220OSONH"].index.intersection(adata_query[adata_query.obs["sampleid"]=="GSM7553445"].obs.index))
#4186

len(pnas.loc[pnas["biosample_id"]=="Hu235OSONH"].index.intersection(adata_query[adata_query.obs["sampleid"]=="GSM7553447"].obs.index))
#2077
len(pnas.index.intersection(adata_query[adata_query.obs["leiden"]=="41"].obs.index))
#0
len(pnas.index.intersection(adata_query[adata_query.obs["leiden"]=="45"].obs.index))
#2
