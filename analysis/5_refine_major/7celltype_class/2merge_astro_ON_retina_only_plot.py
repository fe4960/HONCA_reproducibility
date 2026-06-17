import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astro_ret_ONONH_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_scvi_trg.h5ad")
sc.pl.umap(adata,color="subclass",palette="tab20",save="_Astro_ret_ONONH_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_subclass.png")

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astro_ret_ONONH_mg_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_scvi_trg.h5ad")
dt={"Astro_ONH_GABBR2+":"Astro_ONH1_GABBR2+", "Astro_retina1_PAX5+GABBR2+": "Astro_ONH2_PAX5+GABBR2+", "Astro_retina2_PAX5+ME1+": "Astro_retina1_PAX5+ME1+", "Astro_retina3_NLGN1+": "Astro_retina2_NLGN1+"}
adata.obs["subclass"]=adata.obs["subclass"].replace(dt)

sc.pl.umap(adata,color="subclass",palette="tab20",save="_Astro_ret_ONONH_mg_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_subclass.png", frameon=False)


