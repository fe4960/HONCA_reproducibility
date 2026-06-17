import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import numpy as np
import sys

adata_list=[]

sc.settings.figdir="HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"
dir1="HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"
region=sys.argv[1]
idx=sys.argv[2]
start=1
if "LC" in region:
    start=3

for i in range(start, 6):
    h5ad=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_merge_cutoff20_major_{i}.h5ad"
    adata=sc.read(h5ad)
    csv=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/xenium_spatial_major_{i}.csv"
    spatial=pd.read_csv(csv,index_col=0,header=0)
    adata.obsm["spatial"]=np.array(spatial.loc[adata.obs.index,], dtype=float)
    t=i-2
    bc=pd.read_csv(f"{dir1}/neighb/Selection_{idx}_cells_stats_ONH_{t}_{region}.csv", header=0, comment="#")
    idx1=adata.obs.index.intersection(bc["Cell ID"])
    adata=adata[idx1] 
    txt=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_{i}.txt"
    anno=pd.read_csv(txt, index_col=1,header=0, sep="\t")
    adata.obs["harmony_anno"]=anno.loc[adata.obs.index,"harmony_anno"]
    adata.obs["library_key"]=anno.loc[adata.obs.index,"slide_id"]
    adata_list.append(adata)

adata_full=ad.concat(adata_list)


lab=["harmony_anno" , "library_key"]
for lb in lab:
    adata_full.obs[lb] = adata_full.obs[lb].astype("category")


sq.gr.spatial_neighbors(adata_full, coord_type="generic", delaunay=True, library_key="library_key")

sq.gr.centrality_scores(adata_full, cluster_key="harmony_anno")

# Generate a palette with as many colors as unique categories
sc.pl.set_rcParams_defaults()
palette = sc.pl.palettes.default_20[:adata_full.obs["harmony_anno"].nunique()]

# Assign it
adata_full.uns["harmony_anno_colors"] = palette

sq.pl.centrality_scores(adata_full, cluster_key="harmony_anno", figsize=(16, 5),   save=f"_centrality_score_ONONH_comb_xenium_merge_cutoff20_major_{region}.png" )

sq.gr.co_occurrence(
    adata_full,
    cluster_key="harmony_anno"
)

for cell in adata_full.obs["harmony_anno"].unique():
    sq.pl.co_occurrence(
        adata_full,
        cluster_key="harmony_anno",
        clusters=cell,
        figsize=(10, 10),
        save=f"_{cell}_co_occurrence_ONONH_comb_xenium_merge_cutoff20_major_{region}.png"
    )


sq.gr.nhood_enrichment(adata_full, cluster_key="harmony_anno", library_key="library_key")

adata_full.write(f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_squidpy_{region}.h5ad")


fig, ax = plt.subplots(1, 1, figsize=(7, 7))
sq.pl.nhood_enrichment(
    adata_full,
    cluster_key="harmony_anno",
    library_key="library_key",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax,
    save=f"_nhood_ONONH_comb_xenium_merge_cutoff20_major_{region}.png"
)


#sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])


