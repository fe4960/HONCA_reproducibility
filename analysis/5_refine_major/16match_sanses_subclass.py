import scanpy as sc
import sys
dir="/dfs3b/ruic20_lab/junw42/"
bname=sys.argv[1]
celltype=sys.argv[2]
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"

adata_query1=sc.read_h5ad(f"{outdir}/{bname}.h5ad")

from matplotlib import rcParams
#rcParams["figure.figsize"] = (5,5)
sc.settings.figdir = f"{outdir}/figures"

import pandas as pd
import numpy as np
pnas=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/SCP_UMAPS.csv",header=0)
samplelist=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONH_ON_public_meta.csv'
meta=pd.read_csv(samplelist,header=0,sep=",")
dt=meta[["sample","sampleid"]].set_index("sample").to_dict()["sampleid"]
pnas["bc"]=pnas["NAME"].replace(dt, regex=True)


#adata_query.obs["majorclass_sanes"]="chen"
pnas=pnas.set_index("bc")
idx=pnas.index.intersection(adata_query1.obs.index)

adata_query2=adata_query1[idx]

adata_query2.obs["subclass_sanes"]=pnas.loc[idx,"Type"]

df = adata_query2.obs.groupby(["subclass_sanes", "subclass"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

import matplotlib.pyplot as plt

fig=plt.figure(figsize=(15, 15))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("PMID_37566633")
plt.ylabel("This_study")
fig.savefig(f'{outdir}/figures/{bname}_sanes_current_conf.pdf')


df.to_csv(f"{outdir}/{bname}_sanes_current_conf.csv")

df1 = adata_query2.obs.groupby(["subclass_sanes", "subclass"]).size()
df1.to_csv(f"{outdir}/{bname}_sanes_compare.csv")
