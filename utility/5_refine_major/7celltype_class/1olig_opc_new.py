import scanpy as sc
import anndata as ad
opc=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat.h5ad")

oligo=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat.h5ad")
#oligo=oligo[oligo.obs["leiden"]!="5"]

#oligo.obs["subclass"]=oligo.obs["leiden"]

adata_list=[opc, oligo]
adata=ad.concat(adata_list)
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_rpe_osf_tmcb_ononh.h5ad")
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc.h5ad")
adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_seurat_new.h5ad")
