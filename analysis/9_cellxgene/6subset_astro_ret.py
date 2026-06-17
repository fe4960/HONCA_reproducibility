import scanpy as sc
celltype="Astrocyte"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
adata=sc.read_h5ad(f"{outdir}/{celltype}_subclass_new2.h5ad")
adata=adata[adata.obs["subclass"].isin(["Astro_retina_NLGN1+","Astro_retina","Astro_ONH"])]
adata.write(f"{outdir}/{celltype}_subclass_new2_ret_onh.h5ad" )
