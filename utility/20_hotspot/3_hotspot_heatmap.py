import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib
from scipy.cluster.hierarchy import leaves_list
import matplotlib.pyplot as plt
import sys
plt.rcParams["font.size"] = 30
dir=sys.argv[1]
label=sys.argv[2]
#adata_result = sc.read(
#    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
#)

hs = joblib.load(f"{dir}/{label}_hs.create_modules_years.pkl")
adata = sc.read(f"{dir}/{label}_hs.adata_years.h5ad")

#adata.obsm["X_umap"] = adata_result[adata.obs.index].obsm["X_umap"]

modules = hs.create_modules(min_gene_threshold=160, core_only=True, fdr_threshold=0.01)

#hs.modules = hs.modules.replace({3: 1, 4: 2, 2: 3})
#for x in [
#    "GNB4",
#    "TEAD4",
#    "HAUS2",
#    "VKORC1L1",
#    "LMNB1-DT",
#    "PCGF6",
#    "EIF2AK3",
#    "FAM185A",
#    "AC122719.3",
#    "NEDD4",
#    "LRRC8B",
#    "SPDYA",
#    "AC092910.3",
#    "SV2C",
#    "SOX11",
#    "ABHD17B",
#]:
#    hs.modules[x] = 1
hs.local_correlation_z.to_csv(f"{dir}/{label}_plot_local_correlations.csv")
hs.plot_local_correlations()
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    f"{dir}/{label}_plot_local_correlations.tiff",
    dpi=600,
    transparent=True,
    bbox_inches="tight",
)

ii = leaves_list(hs.linkage)

mod_reordered = hs.modules.iloc[ii]

mod_map = {}
y = np.arange(modules.size)

for x in mod_reordered.unique():
    if x == -1:
        continue

mod_map[x] = y[mod_reordered == x].mean()

mod_reordered.to_csv(
    f"{dir}/{label}_mod_reordered.csv"
)

module_scores = hs.calculate_module_scores()

plt.rcParams["font.size"] = 10


for x in mod_reordered.unique():
    if x == -1:
        continue
    adata.obs[f"Module{x}"] = module_scores.loc[adata.obs.index, x].values
    plt.clf()
    sc.pl.umap(adata, color=f"Module{x}", vmax=6, frameon=False)
    fig = plt.gcf()
    plt.ylabel("")
    plt.xlabel("")
    plt.title(f"Module {x}", fontsize=40)
    fig.set_size_inches(5, 5)
    plt.savefig(
        f"{dir}/{label}_Module{x}_gene_score.tiff",
        dpi=300,
        transparent=True,
        bbox_inches="tight",
    )

    

adata.obs.to_csv(f"{dir}/{label}_Module_gene_score.csv")


###############
#adata_result = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/hotspot_result/PRPC/PRPC_hs.adata_velocity_pseudotime.h5ad")

#sc.pp.normalize_total(adata)
#sc.pp.scale(adata)

for i in mod_reordered.unique():
    if i == -1:
        continue
    top_genes = mod_reordered[mod_reordered==i][:50].index
    scv.pl.heatmap(adata, var_names=top_genes, sortby='age_year',color_map = "RdBu_r")
    fig = plt.gcf()
    fig.set_size_inches(10,4)
    plt.savefig(
         f"{dir}/{label}_gene_module_{i}_heatmap.tiff",
            #        output_file_path + f"{label}_gene_module_{i}_heatmap.tiff",
        bbox_inches="tight",
        transparent=True,
        dpi = 300
    )


