import scanpy as sc
import anndata as ad

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/GSE267301/GSE267301_hvg_10000_epoch_none_seurat_rs_1_scvi_trg.h5ad")
adata_full=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/GSE267301.h5ad")
astro=adata[adata.obs["leiden"].isin(["0","4","7","18"])].obs.index

adata_astro=adata_full[astro].copy()
adata_astro.obs["majorclass"]="Astrocyte"
adata_astro.obs["subclass"]="Astro_Brain"


h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2.h5ad"

adata1=sc.read(h5ad)
adata1.X=adata1.layers["counts"]

adata_full1=ad.concat([adata_astro,adata1])


adata_full1.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge.h5ad")

del adata_astro
del adata1
del adata_full1


astro=adata[adata.obs["leiden"].isin(["1","2","3","12","26"])].obs.index

adata_astro=adata_full[astro].copy()
adata_astro.obs["majorclass"]="Oligodendrocyte"
adata_astro.obs["subclass"]="Oligo_Brain"



h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_new.h5ad"

adata1=sc.read(h5ad)
adata1.X=adata1.layers["counts"]

adata_full1=ad.concat([adata_astro,adata1])


adata_full1.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge.h5ad")

