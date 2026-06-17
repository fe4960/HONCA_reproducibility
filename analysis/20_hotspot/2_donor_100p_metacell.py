#import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib
import sys
from sklearn.cluster import KMeans

h5ad=sys.argv[1]
label=sys.argv[2]
t=sys.argv[3]
output_path = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/"
adata = sc.read(h5ad)

#adata = adata[adata.obs.majorclass.isin(["PRPC"])]

adata.obs["age_year"] = adata.obs["age_year"].astype(float)
adata = adata[adata.obs["age_year"] > 0]



#sc.pp.calculate_qc_metrics(adata, inplace=True)
#adata = adata[:, adata.var.mean_counts > 0]

del adata.X


adata.X=adata.layers["counts"].copy()

for dn in adata.obs[t].unique()
    adata_t=adata[adata.obs[t]==dn]
    
    n_cells=adata_t.n_obs
    n_per_meta=50
    n_meta = n_cells // n_per_meta

    np.random.seed(7)
    shuffled_indices = np.random.permutation(n_cells)

    meta_exprs = []
    meta_obs = []

    for i in range(n_meta):
        idx = shuffled_indices[i*n_per_meta:(i+1)*n_per_meta]
        avg_expr = adata.X[idx].mean(axis=0)
        meta_exprs.append(avg_expr)
        meta_obs.append(f"meta_cell_{i}")

# Step 4: Create a new AnnData object
    meta_exprs = np.array(meta_exprs)
    meta_adata = sc.AnnData(meta_exprs)
    meta_adata.var = adata.var.copy()
    meta_adata.obs_names = meta_obs

# Step 5: Save or analyze
    meta_adata.write("meta_cells_100.h5ad")

    print(f"Generated {meta_adata.n_obs} meta cells from {n_cells} original cells.")


#sc.pp.highly_variable_genes(
#        adata, flavor="seurat_v3", n_top_genes=10000, subset=True
#    )

# Step 4: Run k-means clustering
n_clusters = 100  # You can change this based on your desired number of meta cells
kmeans = KMeans(n_clusters=n_clusters, random_state=0)
adata.obs['kmeans'] = kmeans.fit_predict(adata.obsm['X_scVI']).astype(str)


meta_cells = adata.to_df().groupby([adata.obs['kmeans'],adata.obs['donor']]).sum()

meta_adata = sc.AnnData(meta_cells)
meta_adata.obs['kmeans'] = meta_cells.index



#sc.pp.subsample(adata, n_obs=10000, random_state=0)

#adata.layers["counts"] = adata.X

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]
adata

###adata.obsm["latent_time"] = np.asarray([[x] for x in adata.obs["latent_time"]])
adata.obsm["age_year"] = np.asarray([[x] for x in adata.obs["age_year"]])

hs = hotspot.Hotspot(
    adata,
    model="danb",
    layer_key="counts",
    latent_obsm_key="age_year",
    umi_counts_obs_key="total_counts",
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=12)

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index  # Select genes
local_correlations = hs.compute_local_correlations(
    hs_genes, jobs=12
)  # jobs for parallelization

modules = hs.create_modules(min_gene_threshold=160, core_only=True, fdr_threshold=0.05)
module_scores = hs.calculate_module_scores()

local_correlations.to_csv(
    output_path + "/" + label + "_hs_local_correlations_years.csv"
)
module_scores.to_csv(
    output_path + "/" + label + "_hs_module_scores_years.csv"
)
joblib.dump(
    hs,
    output_path + "/" + label + "_hs.create_modules_years.pkl",
)
adata.write(output_path +  "/" + label +  "_hs.adata_years.h5ad"
)
