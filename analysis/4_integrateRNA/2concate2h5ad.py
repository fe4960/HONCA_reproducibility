import scanpy as sc
import numpy as np
import sys
import pandas as pd
import anndata as ad

samplelist=sys.argv[1]
bname=sys.argv[2]
outdir=sys.argv[3]

sl=pd.read_csv(samplelist,sep="\t",header=0)

adata_list=[]

#adata_full=ad.AnnData()

#i=0
#adata=sc.read_h5ad(sl.loc[i,sl.columns[-1]])



for i in sl.index :
    adata=sc.read_h5ad(sl.loc[i,sl.columns[-1]])
    for c in sl.columns :
        if c != "h5ad":
            adata.obs[c]=sl.loc[i,c]
    adata.obs.index=sl.loc[i,sl.columns[0]]+"_"+adata.obs.index
#adata_full = adata_full.concatenate(adata, index_unique=None)
    adata_list.append(adata)


adata_full=ad.concat(adata_list)

adata_full.write(f"{outdir}/{bname}.h5ad")


