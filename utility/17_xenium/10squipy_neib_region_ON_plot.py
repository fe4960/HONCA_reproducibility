import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import numpy as np
import sys


sc.settings.figdir="HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"
dir1="HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"
#region=sys.argv[1]
res=["before_LC", "after_LC", "ON", "mid"]

for region in res:
    adata_full=sc.read(f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_squidpy_{region}_clean.h5ad")
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    sq.pl.nhood_enrichment(
        adata_full,
        cluster_key="harmony_anno",
        library_key="library_key",
        figsize=(5, 5),
        title="Neighborhood enrichment adata",
        ax=ax,
        save=f"_nhood_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean_st.png",
        vmax=54.35,
        vmin=-28.41
        )
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    sq.pl.nhood_enrichment(
        adata_full,
        cluster_key="harmony_anno",
        library_key="library_key",
        figsize=(5, 5),
        title="Neighborhood enrichment adata",
        ax=ax,
        save=f"_nhood_ONONH_comb_xenium_merge_cutoff20_major_{region}_clean_cen.png",
        vmax=30,
        vmin=-30
    )
