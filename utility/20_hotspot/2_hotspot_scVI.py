#import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib
import sys

h5ad=sys.argv[1]
label=sys.argv[2]
hvg=int(sys.argv[3])
output_path = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/"
adata = sc.read(h5ad)


#adata.obs["age_year"] = adata.obs["age_year"].astype(float)
#adata = adata[adata.obs["age_year"] > 0]


sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]

del adata.X


adata.X=adata.layers["counts"].copy()

sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=hvg, subset=True
    )


sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:, adata.var.mean_counts > 0]
adata


adata.obs["total_counts"] = np.asarray(adata.X.sum(1)).ravel()

###adata.obsm["latent_time"] = np.asarray([[x] for x in adata.obs["latent_time"]])
#adata.obsm["age_year"] = np.asarray([[x] for x in adata.obs["age_year"]])

hs = hotspot.Hotspot(
    adata,
    model="danb",
    layer_key="counts",
#    latent_obsm_key="age_year",
    latent_obsm_key="X_scVI",
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
    output_path + "/" + label + f"_hs_local_correlations_scvi_hvg_{hvg}.csv"
)
module_scores.to_csv(
    output_path + "/" + label + f"_hs_module_scores_scvi_hvg_{hvg}.csv"
)
joblib.dump(
    hs,
    output_path + "/" + label + f"_hs.create_modules_scvi_hvg_{hvg}.pkl",
)
adata.write(output_path +  "/" + label +  f"_hs.adata_scvi_hvg_{hvg}.h5ad"
)
