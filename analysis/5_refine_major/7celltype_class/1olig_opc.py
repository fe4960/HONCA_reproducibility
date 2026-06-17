import scanpy as sc
import anndata as ad
opc=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass.h5ad")

#oligo=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_res_1_clean_res_1.0.h5ad")
oligo=sc.read("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_scvi_trg.h5ad")
#oligo=oligo[-oligo.obs["leiden"].isin(["15","16","4"])]
oligo=oligo[oligo.obs["leiden"]!="5"]

oligo.obs["subclass"]=oligo.obs["leiden"]

adata_list=[opc, oligo]
adata=ad.concat(adata_list)
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_rpe_osf_tmcb_ononh.h5ad")
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc.h5ad")
adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_new.h5ad")
