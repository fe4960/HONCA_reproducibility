#import scvelo as scv
#import system 
import scanpy as sc
import anndata
import numpy as np
import sys
#adata=scv.read("/storage/chenlab/Users/junwang/human_meta/data/ref/_major_b_rgc_sample_major_scvi_Cluster.h5ad")
fn=sys.argv[1]
cluster_key=sys.argv[2]
target_cells = int(sys.argv[3])
out=sys.argv[4]

label=sys.argv[5]
celltype=sys.argv[6]

adata=sc.read(f'{fn}.h5ad')

if celltype != "major":
    adata=adata[adata.obs[label]==celltype].copy()

adata.obs[cluster_key]=adata.obs[cluster_key].astype('category')
#adatas = [adata[adata.obs[cluster_key].isin(clust)] for clust in adata.obs[cluster_key].cat.categories]
adatas=[adata[adata.obs[cluster_key] == clust] for clust in adata.obs[cluster_key].cat.categories]

for dat in adatas:
    if dat.n_obs > target_cells:
         sc.pp.subsample(dat, n_obs=target_cells)

adata_downsampled=anndata.concat(adatas[0:])

adata_downsampled.write(f'{out}.h5ad')

adata_downsampled.obs.to_csv(f'{out}.obs.txt.gz',sep="\t")

