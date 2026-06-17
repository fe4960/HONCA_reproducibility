import scanpy as sc

adata=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]

del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_clusters.h5ad")

adata=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]
del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling_clusters.h5ad")



adata=sc.read("HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]
del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_clusters.h5ad")



adata=sc.read("HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]
del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location_clusters.h5ad")

