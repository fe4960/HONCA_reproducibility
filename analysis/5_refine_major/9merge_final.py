import scanpy as sc
import anndata as ad
ct=["Oligodendrocyte","Fibroblast","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","HC","Cone","NK_T_cell","RGC","Schwann_cell"]

ad_list=[]

for c in ct:
    adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{c}/clean/{c}_subclass.h5ad')
    ad_list.append(adata)

adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new.h5ad')
ad_list.append(adata)

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Melanocyte/clean/Melanocyte_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
adata.obs["subclass"]="Melanocyte"
ad_list.append(adata)

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Rod/Rod_hvg2k_epochnone_seurat_v3_zinb_scvi_trg.h5ad")
adata.obs["subclass"]="Rod"
ad_list.append(adata)


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/B_cell/B_cell_hvg2k_epochnone_seurat_v3_zinb_scvi_trg.h5ad")
adata.obs["subclass"]="Adipocyte"
ad_list.append(adata)

adata_full=ad.concat(ad_list)
adata_full.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/ON_ONH_clean.h5ad")
