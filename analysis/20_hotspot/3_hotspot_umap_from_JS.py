import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import mplscience
dir=sys.argv[1]
label=sys.argv[2]
t=sys.argv[3]
hvg=sys.argv[4]

######hs = joblib.load(f"{dir}/{label}_hs.create_modules_years.pkl")
#######adata = sc.read(f"{dir}/{label}_hs.adata_years.h5ad")

sc.settings.figdir=dir

hs = joblib.load(f"{dir}/{label}_hs.create_modules_{t}_{hvg}.pkl")
adata = sc.read(f"{dir}/{label}_hs.adata_{t}_{hvg}.h5ad")


 
#adata.obs.age=pd.to_numeric(adata.obs[t], errors='coerce')
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

#for min_gene in [100,120, 160,200, 250, 300]:
for min_gene in [200]:
    print(min_gene)
    modules = hs.create_modules(min_gene_threshold=min_gene, core_only=True, fdr_threshold=0.05)
    module_scores = hs.calculate_module_scores()
    module_scores.to_csv(f"{dir}/{label}_hs_gene_module_scores_min_gene_{min_gene}_{hvg}_{t}.csv")
    #plot loca correlation
    hs.plot_local_correlations()
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    plt.savefig(
        f"{dir}/{label}_min_gene_{min_gene}_plot_local_correlations_{hvg}_{t}.tiff",
        dpi=300,
        transparent=True,
        bbox_inches="tight",
    )
    #output gene module assignment
    results = hs.results.join(hs.modules)
    results.to_csv(f"{dir}/{label}_hs_gene_module_assignment_min_gene_{min_gene}_{hvg}_{t}.csv")

#    ii = leaves_list(hs.linkage)
    mod_reordered = hs.modules

    mod_map = {}
    y = np.arange(modules.size)

    for x in mod_reordered.unique():
        if x == -1:
            continue

        mod_map[x] = y[mod_reordered == x].mean()

    mod_reordered.to_csv(
    f"{dir}/{label}_mod_reordered_{min_gene}_{hvg}_{t}.csv"
    )

    module_scores = hs.calculate_module_scores()

    plt.rcParams["font.size"] = 10

    module_cols = []
    for x in hs.modules.unique():
        if x == -1 :
            continue
        adata.obs[f"Module{x}"] = module_scores.loc[adata.obs.index, x].values
        module_cols.append(f"Module{x}")

    #module_cols = ["Module1", "Module2", "Module3", "Module4", "Module5", "Module6", "Module7"]
#    module_cols = ["Module1", "Module2", "Module3", "Module4", "Module5"]
    module_cols = ["Module1", "Module2", "Module3", "Module4", "Module5", "Module6"]
    with mplscience.style_context():
        sc.pl.umap(adata, color=module_cols, frameon=False, vmin=-1, vmax=1, save=f"{label}_min_gene_{min_gene}_{t}_{hvg}.png", ncols=3)

#        sc.pl.umap(adata, color=module_cols, frameon=False, vmin=-1, vmax=1, save=f"{label}_min_gene_{min_gene}_{t}_{hvg}.png")

