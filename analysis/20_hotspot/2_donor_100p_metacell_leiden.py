import scanpy as sc
import numpy as np
import pandas as pd
import sys
h5ad=sys.argv[1]
label=sys.argv[2]
adata = sc.read_h5ad(f"{h5ad}.h5ad")

del adata.X

adata.X=adata.layers["counts"].copy()

sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_scVI")
sc.tl.leiden(adata, resolution=adata.n_obs / 1000)
adata.obs['meta_cluster']=adata.obs['leiden']

# Step 4: Aggregate expression by cluster
meta_exprs = []
meta_obs = []
meta_dn = []
meta_age = []
for dn in adata.obs[label].unique():
    adata1=adata[adata.obs[label]==dn]
    age = adata1.obs["age_year"].unique()[0]
    for cluster_id in sorted(adata1.obs['meta_cluster'].unique(), key=int):
        idx = adata1.obs['meta_cluster'] == cluster_id
        X=adata1[idx].X
        avg_expr = np.ravel(X.sum(axis=0))
        meta_exprs.append(avg_expr)
        meta_obs.append(f"metacell_{cluster_id}_{dn}")
        meta_dn.append(dn)
        meta_age.append(age)

# Step 5: Create new AnnData for meta-cells
meta_exprs = np.array(meta_exprs)
meta_adata = sc.AnnData(meta_exprs)
meta_adata.var = adata.var.copy()
meta_adata.obs_names = meta_obs
meta_adata.obs['meta_cluster'] = meta_adata.obs_names
meta_adata.obs['age']=meta_age
meta_adata.obs['donor']=meta_dn

meta_adata.write(f"{h5ad}_leiden_{label}_metacell.h5ad")

