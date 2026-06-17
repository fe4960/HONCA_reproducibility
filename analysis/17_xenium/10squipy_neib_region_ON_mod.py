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
#idx3=sys.argv[3]
start=1
#if "LC" in region:
#    start=3


if "ONL" in region:
    for i in range(1, 3):
        h5ad=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_merge_cutoff20_major_{i}.h5ad"
        adata=sc.read(h5ad)
        csv=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/xenium_spatial_major_{i}.csv"
        spatial=pd.read_csv(csv,index_col=0,header=0)
        adata.obsm["spatial"]=np.array(spatial.loc[adata.obs.index,], dtype=float)
        adata=adata[adata.obs["nFeature_RNA"]>=5].copy()
#    t=i-2
        bc=pd.read_csv(f"{dir1}/neighb/Selection_{idx}_cells_stats_PP_{i}_{region}.csv", header=0, comment="#")
        idx1=adata.obs.index.intersection(bc["Cell ID"])
        adata=adata[idx1] 
        txt=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_{i}.txt"
        anno=pd.read_csv(txt, index_col=1,header=0, sep="\t")
        adata.obs["harmony_anno"]=anno.loc[adata.obs.index,"harmony_anno"]

        adata.obs["library_key"]=anno.loc[adata.obs.index,"slide_id"]
        unique_ids = anno["slide_id"].unique()
        slide_id = unique_ids[0] if len(unique_ids) > 0 else None
        csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/astro_{slide_id}_subclass1_RCTD_results_nonneuron.txt"
#        csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/astro_{slide_id}_RCTD_results_nonneuron.txt"

        oligo=pd.read_csv(csv,index_col=0,header=0, sep="\t")
        idx2=oligo.index.intersection(adata.obs.index)
        adata.obs.loc[idx2,"harmony_anno"]=oligo.loc[idx2, "first_type"]
        del_idx=oligo.loc[(oligo["spot_class"]=="reject")| (oligo["spot_class"]=="doublet_uncertain")].index.copy()
        csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_{slide_id}_subclass_RCTD_results_nonneuron01.txt"

#        csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_{slide_id}_subclass_RCTD_results_nonneuron015.txt"
        oligo1=pd.read_csv(csv,index_col=0,header=0, sep="\t")
        idx2=oligo1.index.intersection(adata.obs.index)
        adata.obs.loc[idx2,"harmony_anno"]=oligo1.loc[idx2, "first_type"]
        del_idx1=oligo1.loc[(oligo1["spot_class"]=="reject")| (oligo1["spot_class"]=="doublet_uncertain")].index.copy()
        adata=adata[~adata.obs.index.isin(del_idx) & ~adata.obs.index.isin(del_idx1) ].copy()
        adata_list.append(adata)



for i in range(3, 6):
    h5ad=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_merge_cutoff20_major_{i}.h5ad"
    adata=sc.read(h5ad)
    csv=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/xenium_spatial_major_{i}.csv"
    spatial=pd.read_csv(csv,index_col=0,header=0)
    adata.obsm["spatial"]=np.array(spatial.loc[adata.obs.index,], dtype=float)
    adata=adata[adata.obs["nFeature_RNA"]>=5].copy()
    t=i-2
    bc=pd.read_csv(f"{dir1}/neighb/Selection_{idx}_cells_stats_ONH_{t}_{region}.csv", header=0, comment="#")
    idx1=adata.obs.index.intersection(bc["Cell ID"])
    adata=adata[idx1] 
    txt=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_{i}.txt"
    anno=pd.read_csv(txt, index_col=1,header=0, sep="\t")
    adata.obs["harmony_anno"]=anno.loc[adata.obs.index,"harmony_anno"]
    adata.obs["library_key"]=anno.loc[adata.obs.index,"slide_id"]
    unique_ids = anno["slide_id"].unique()
    slide_id = unique_ids[0] if len(unique_ids) > 0 else None
    csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/astro_{slide_id}_subclass1_RCTD_results_nonneuron.txt"
    oligo=pd.read_csv(csv,index_col=0,header=0, sep="\t")
    idx2=oligo.index.intersection(adata.obs.index)
    adata.obs.loc[idx2,"harmony_anno"]=oligo.loc[idx2, "first_type"]
    del_idx=oligo.loc[oligo["spot_class"]!="singlet"].index.copy()
    csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_{slide_id}_subclass_RCTD_results_nonneuron01.txt"

#    csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_{slide_id}_subclass_RCTD_results_nonneuron015.txt"
    oligo1=pd.read_csv(csv,index_col=0,header=0, sep="\t")
    idx2=oligo1.index.intersection(adata.obs.index)
    adata.obs.loc[idx2,"harmony_anno"]=oligo1.loc[idx2, "first_type"]
    del_idx1=oligo1.loc[oligo1["spot_class"]!="singlet"].index.copy()
    adata=adata[~adata.obs.index.isin(del_idx) & ~adata.obs.index.isin(del_idx1) ].copy()

    adata_list.append(adata)

adata_full=ad.concat(adata_list)

if region == "before_LC" :
    adata_full=adata_full[~adata_full.obs["harmony_anno"].isin(["Rod","Cone","BC","Schwann_cell","RGC","Melanocyte","MG", "HC", "Astrocyte", "Oligodendrocyte", "astro_ON", "OPC", "astro_ONHON"])] # "OPC", "Oligodendrocyte", "HC"])]
elif region == "after_LC" :
    adata_full=adata_full[~adata_full.obs["harmony_anno"].isin(["Rod","Cone","BC","Schwann_cell","RGC","Melanocyte","MG", "HC", "Astrocyte", "Oligodendrocyte", "RPE"])]
elif region in [ "mid", "ON"]  :
    adata_full=adata_full[~adata_full.obs["harmony_anno"].isin(["Rod","Cone","BC","Schwann_cell","RGC","Melanocyte","MG", "HC", "Astrocyte", "Oligodendrocyte", "astro_retina", "astro_ONH1", "astro_ONH2", "RPE"])]


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

sq.pl.centrality_scores(adata_full, cluster_key="harmony_anno", figsize=(16, 5),   save=f"_centrality_score_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean_refine.png" )

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
        save=f"_{cell}_co_occurrence_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean_refine.png"
    )


sq.gr.nhood_enrichment(adata_full, cluster_key="harmony_anno", library_key="library_key")

adata_full.write(f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_squidpy_{region}_clean_refine.h5ad")


fig, ax = plt.subplots(1, 1, figsize=(7, 7))
sq.pl.nhood_enrichment(
    adata_full,
    cluster_key="harmony_anno",
    library_key="library_key",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax,
    #cmap="seismic",
   # save=f"_nhood_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean.png",
    save=f"_nhood_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean_st_refine.pdf",
#    save=f"_nhood_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean_cen.png",

#    vmax=54.35,
#    vmin=-28.41
#     vmin=30,

)


#sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])


