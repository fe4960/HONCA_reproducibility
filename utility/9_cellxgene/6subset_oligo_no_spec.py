import scanpy as sc
celltype="Oligodendrocyte"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
adata=sc.read_h5ad(f"{outdir}/{celltype}_subclass_new.h5ad")
adata=adata[-adata.obs["subclass"].isin(["OLIGO_non_specific"])]
adata.write(f"{outdir}/{celltype}_subclass_new_rm_non_spec.h5ad" )
