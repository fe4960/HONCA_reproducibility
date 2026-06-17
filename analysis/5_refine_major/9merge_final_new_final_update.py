import scanpy as sc
import anndata as ad
import pandas as pd
#ct=["Oligodendrocyte","Fibroblast","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","HC","Cone","NK_T_cell","RGC","Schwann_cell"]

ct=["BC","AC","RPE", "HC","Cone","NK_T_cell"]

#ct1=["Endothelial_cell", "Macrophage", "MG", "Microglia", "RGC", "Schwann_cell"]

ct1=[ "Macrophage", "MG", "Microglia", "RGC", "Schwann_cell"]


ad_list=[]

for c in ct:
    adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{c}/clean/{c}_subclass.h5ad')
    ad_list.append(adata)

for c in ct1:
    adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{c}/clean/{c}_subclass_sb.h5ad')
    ad_list.append(adata)



adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad')
ad_list.append(adata)

adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling.h5ad')
ad_list.append(adata)


adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location.h5ad')

#adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new.h5ad')
ad_list.append(adata)

adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_subclass_sb_seurat.h5ad')
ad_list.append(adata)

adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/Mural_cell_subclass_rmRPE.h5ad')
ad_list.append(adata)


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Melanocyte/clean/Melanocyte_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
adata.obs["subclass"]="Melanocyte"
ad_list.append(adata)

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Rod/Rod_hvg2k_epochnone_seurat_v3_zinb_scvi_trg.h5ad")
adata.obs["subclass"]="Rod"
ad_list.append(adata)


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/B_cell/B_cell_hvg2k_epochnone_seurat_v3_zinb_scvi_trg.h5ad")
adata.obs["subclass"]="Adipocyte"
adata.obs["majorclass"]="Adipocyte"
adata.obs["majorclass1"]="Adipocyte"

ad_list.append(adata)

#adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_scvi_trg.h5ad')

#obs=pd.read_csv(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_res_0.3_clean_rm4_15_res_0.3.obs.gz', header=0, index_col=0)

#adata.obs["subclass"]=obs.loc[adata.obs.index,"leiden"].astype(str)

adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat.h5ad")

ad_list.append(adata)


adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]

adata_full=ad.concat(ad_list)
adata_full.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/ON_ONH_clean_new_update_final.h5ad")
adata_full.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/ON_ONH_clean_new_update_final.obs.gz")
adata_full.var.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/ON_ONH_clean_new_update_final.var.gz")
