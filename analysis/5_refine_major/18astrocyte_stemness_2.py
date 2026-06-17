import scanpy as sc

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/figures"

#mk=["SOX2", "NES", "PROM1", "SOX9", "PAX6", "FABP7", "ITGA6",  "NFIA", "NFIX", "KIF2C", "TOP2A", "SLC1A3", "GFAP", "VIM", "RBPJ", "GJB6", "EGFR", "ASCL1", "MKI67", "DCX", "NAV3"]
#adata=sc.read()
#adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad')
#adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_rmRPE.h5ad')
#adata.obs.to_csv(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_rmRPE.obs.gz')
#sc.pl.dotplot(adata,mk, groupby="subclass", save="_Astrocyte_subclass_new5_clean_stemness.png", use_raw=False)

celltype="Oligodendrocyte"
indir="/dfs3b/ruic20_lab/junw42"
adata=sc.read(f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/oligo_opc_seurat_new_scvi_cycling.h5ad")


#sc.pl.umap(adata, color=mk, save="_Astrocyte_subclass_new5_clean_stemness_umap.png")


#NSC_quiescent=["GFAP", "SOX2", "NES", "VIM", "FABP7", "HES1", "HES5", "ID1", "ID2", "ID3", "PROM1", "AQP4", "SLC1A3", "ALDH1L1", "GLI3"]

#sc.pl.dotplot(adata, NSC_quiescent, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NSC_quiescent_umap.png")

#NSC_activated=[ "SOX2", "NES", "VIM", "EGFR", "ASCL1", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM6",   "HMGB2", "TUBB", "HES1", "DLL1", "CCND1"]

#sc.pl.dotplot(adata, NSC_activated, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NSC_activated_umap.png")

#NPC_early=["ASCL1", "DLL1", "EGFR", "SOX4", "SOX11", "DCX", "ELAVL4", "STMN2", "TUBB3", "NEUROD1", "NEUROD2", "HES6", "BTG2", "MKI67", "TOP2A"]

#sc.pl.dotplot(adata, NPC_early, groupby="subclass", save="_Astrocyte_subclass_new5_clean_NPC_early_umap.png")



#Astrocyte_homeostatic=["AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1",  "FGFR3", "ATP1A2", "GLUL", "APOE", "SOX9"]
#dt={"Astro_homo": ["AQP4", "SLC1A2", "SLC1A3", "GJA1",  "FGFR3", "APOE", "SOX9"], "Activated and progenitor-like": ["GFAP", "SOX2", "VIM", "ID1", "SOX11", "HES6", "NES", "TUBB"]}
#mk=["AQP4", "APOE", "SOX9", "GFAP", "SOX2", "VIM", "GJA1",  "ID1", "SOX11", "HES6", "NES", "TUBB"]
#sc.pl.dotplot(adata, mk, groupby="subclass", save="_Astrocyte_subclass_new5_clean_stemness_important.pdf")

level=["OPC", "OPC_GLIS3+", "OPC_Cycling", "COP_NFOL", "COP_LAMA2+", "MFOL", "OLIGO2_LRRC7+", "OLIGO1", "OLIGO1_SVEP1+", "OLIGO1_RBFOX1+"]

adata.obs["subclass"] = (
    adata.obs["subclass"]
    .astype("category")
    .cat.set_categories(level, ordered=True)
)

mk=["PIEZO2", "TRPC1", "PKD2", "TRPM7"]
sc.pl.dotplot(adata, mk, groupby="subclass", save="_Oligo_opc_subclass_seurat_rmRPE_mechano.pdf", swap_axes=True)

level=["OPC", ""]

mk=["PIEZO2"]
sc.pl.dotplot(adata, mk, groupby="subclass", save="_Oligo_opc_subclass_seurat_rmRPE_mechano1.pdf", categories_order=level, swap_axes=True)

