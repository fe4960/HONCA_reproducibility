import scanpy as sc




adata=sc.read("HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_subclass_rmRPE.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]

del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_subclass_rmRPE_clusters.h5ad")

adata=sc.read("HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/Mural_cell_subclass_rmRPE.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]
del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/Mural_cell_subclass_rmRPE_clusters.h5ad")



adata=sc.read("HCA_ON/data/5_refine_major/scvi/Microglia/clean/Microglia_subclass_sb_clean.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]
del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Microglia/clean/Microglia_subclass_sb_clean_clusters.h5ad")



adata=sc.read("HCA_ON/data/5_refine_major/scvi/Macrophage/clean/Macrophage_subclass_sb_clean_rmRPE.h5ad")
adata.obs['clusters']=adata.obs['subclass']
adata.obs['Sample']=adata.obs['sampleid']
adata.X=adata.layers["counts"]
del adata.layers["counts"]

adata.write("HCA_ON/data/5_refine_major/scvi/Macrophage/clean/Macrophage_subclass_sb_clean_rmRPE_clusters.h5ad")

