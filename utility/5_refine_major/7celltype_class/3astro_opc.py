import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astro_opc_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_scvi_trg.h5ad")
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures/"
sc.pl.umap(adata,color="subclass",save="Astro_opc_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_subclass.png")

sc.pl.umap(adata,color="subclass",groups=["OPC_GLIS3+"],save="Astro_opc_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_subclass_GLIS3.png")
