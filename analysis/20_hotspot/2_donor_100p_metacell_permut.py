#import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
#import hotspot
import joblib
import sys
from sklearn.cluster import KMeans

h5ad=sys.argv[1]
label=sys.argv[2]
t=sys.argv[3]
n=sys.argv[4]
output_path = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/"
adata = sc.read(h5ad)

#adata = adata[adata.obs.majorclass.isin(["PRPC"])]

adata.obs["age_year"] = adata.obs["age_year"].astype(float)
adata = adata[adata.obs["age_year"] > 0]



#sc.pp.calculate_qc_metrics(adata, inplace=True)
#adata = adata[:, adata.var.mean_counts > 0]

del adata.X


adata.X=adata.layers["counts"].copy()

meta_exprs = []
meta_obs = []
meta_dn = []
meta_age = []

for dn in adata.obs[t].unique():
    adata_t=adata[adata.obs[t]==dn]
    age=adata_t.obs["age_year"].astype(int).unique()[0]
    n_cells=adata_t.n_obs
    n_per_meta=int(n)
    n_meta = n_cells // n_per_meta

    np.random.seed(7)
    shuffled_indices = np.random.permutation(n_cells)


    for i in range(n_meta):
        idx = shuffled_indices[i*n_per_meta:(i+1)*n_per_meta]
#        avg_expr = adata[idx].X.sum(axis=0)
        X=adata_t[idx].X
        avg_expr = np.ravel(X.sum(axis=0))
        meta_exprs.append(avg_expr)
        meta_obs.append(f"meta_cell_{i}_{dn}")
        meta_dn.append(f"{dn}")
        meta_age.append(age)

# Step 4: Create a new AnnData object
meta_exprs = np.array(meta_exprs)
meta_adata = sc.AnnData(meta_exprs)
meta_adata.var = adata.var.copy()
meta_adata.obs_names = meta_obs
meta_adata.obs["donor"]=meta_dn
meta_adata.obs["age"]=meta_age

# Step 5: Save or analyze
#meta_adata.write("meta_cells_100.h5ad")

meta_adata.write(f"{h5ad}_leiden_{label}_permute_metacell_{n}.h5ad")

#print(f"Generated {meta_adata.n_obs} meta cells from {n_cells} original cells.")


