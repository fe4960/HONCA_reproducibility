import scanpy as sc

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures"

mk=["SOX2", "NES", "PROM1", "SOX9", "PAX6", "FABP7", "ITGA6",  "NFIA", "NFIX", "KIF2C", "TOP2A", "SLC1A3", "GFAP", "VIM", "RBPJ", "GJB6", "EGFR", "ASCL1", "MKI67", "DCX", "NAV3"]
#adata=sc.read()
adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad')

sc.pl.dotplot(adata,mk, groupby="subclass", save="_Astrocyte_subclass_new5_clean_stemness.png", use_raw=False)

sc.pl.umap(adata, color=mk, save="_Astrocyte_subclass_new5_clean_stemness_umap.png")


NSC_quiescent=["GFAP", "SOX2", "NES", "VIM", "FABP7", "HES1", "HES5", "ID1", "ID2", "ID3", "PROM1", "AQP4", "SLC1A3", "ALDH1L1", "GLI3"]

sc.pl.dotplot(adata, NSC_quiescent, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NSC_quiescent_umap.png")

NSC_activated=[ "SOX2", "NES", "VIM", "EGFR", "ASCL1", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM6",   "HMGB2", "TUBB", "HES1", "DLL1", "CCND1"]

sc.pl.dotplot(adata, NSC_activated, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NSC_activated_umap.png")

NPC_early=["ASCL1", "DLL1", "EGFR", "SOX4", "SOX11", "DCX", "ELAVL4", "STMN2", "TUBB3", "NEUROD1", "NEUROD2", "HES6", "BTG2", "MKI67", "TOP2A"]

sc.pl.dotplot(adata, NPC_early, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NPC_early_umap.png")



Astrocyte_homeostatic=["AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1",  "FGFR3", "ATP1A2", "GLUL", "APOE", "SOX9"]
sc.pl.dotplot(adata, Astrocyte_homeostatic, groupby="subclass", save="_Astrocyte_subclass_new5_clean_Astro_homo_umap.png")




