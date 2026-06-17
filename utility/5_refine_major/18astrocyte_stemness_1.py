import scanpy as sc

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures"

#mk=["SOX2", "NES", "PROM1", "SOX9", "PAX6", "FABP7", "ITGA6",  "NFIA", "NFIX", "KIF2C", "TOP2A", "SLC1A3", "GFAP", "VIM", "RBPJ", "GJB6", "EGFR", "ASCL1", "MKI67", "DCX", "NAV3"]
#adata=sc.read()
adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad')

#sc.pl.dotplot(adata,mk, groupby="subclass", save="_Astrocyte_subclass_new5_clean_stemness.png", use_raw=False)

#sc.pl.umap(adata, color=mk, save="_Astrocyte_subclass_new5_clean_stemness_umap.png")


#NSC_quiescent=["GFAP", "SOX2", "NES", "VIM", "FABP7", "HES1", "HES5", "ID1", "ID2", "ID3", "PROM1", "AQP4", "SLC1A3", "ALDH1L1", "GLI3"]

#sc.pl.dotplot(adata, NSC_quiescent, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NSC_quiescent_umap.png")

#NSC_activated=[ "SOX2", "NES", "VIM", "EGFR", "ASCL1", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM6",   "HMGB2", "TUBB", "HES1", "DLL1", "CCND1"]

#sc.pl.dotplot(adata, NSC_activated, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NSC_activated_umap.png")

#NPC_early=["ASCL1", "DLL1", "EGFR", "SOX4", "SOX11", "DCX", "ELAVL4", "STMN2", "TUBB3", "NEUROD1", "NEUROD2", "HES6", "BTG2", "MKI67", "TOP2A"]

#sc.pl.dotplot(adata, NPC_early, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NPC_early_umap.png")




level=[  "Astro_retina1_PAX5+ME1+", "Astro_retina2_NLGN1+", "Astro_ONH1_GABBR2+SV2B+", "Astro_ONH2_PAX5+GABBR2+", "Astro_ONHON1_SLC4A11+MARCH1+", "Astro_ONHON2_CST3+APOE+",  "Astro_ON1_WNK2+TSHZ2+","Astro_ON2_ACTN1+SERPINA3+", "Astro_ON3_NR4A3+FOS+", "Astro_ON4_AFF3+DMGDH+", "Astro_ON5_DPP10+"]


adata.obs["subclass"] = (
    adata.obs["subclass"]
    .astype("category")
    .cat.set_categories(level, ordered=True)
)

#Astrocyte_homeostatic=["AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1",  "FGFR3", "ATP1A2", "GLUL", "APOE", "SOX9"]
#dt={"Astro_homo": ["AQP4", "SLC1A2", "SLC1A3", "GJA1",  "FGFR3", "APOE", "SOX9"], "Activated and progenitor-like": ["GFAP", "SOX2", "VIM", "ID1", "SOX11", "HES6", "NES", "TUBB"]}
mk=["AQP4", "APOE", "SOX9", "GFAP", "SOX2", "VIM", "GJA1",  "ID1", "SOX11", "HES6", "NES", "TUBB"]
sc.pl.dotplot(adata, mk, groupby="subclass", save="_Astrocyte_subclass_new5_clean_stemness_important.pdf", categories_order=level)

mk=["PIEZO2", "PIEZO1"]

sc.pl.dotplot(adata, mk, groupby="subclass", save="_Astrocyte_subclass_new5_clean_mechano.pdf", swap_axes=True, categories_order=level)



