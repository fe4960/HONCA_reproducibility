import scanpy as sc
import pandas as pd
adata=sc.read("/dfs3b/ruic20_lab/jianms1/Single_cell/Paneye/Endo/Endo_snRNA_all/annotation/snRNA_metaEndo_paneye_leiden_1_0.9_0.05_annotated_normexp.h5ad")

adata_on=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_0.6_clean_rm_sb_none_scvi_trg.h5ad")

idx=adata_on.obs.index.intersection(adata.obs.index)
adata_on.obs.loc[idx,"anno"]=adata.obs.loc[idx,'subclass_meta_endo']

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/figures"

sc.pl.umap(adata_on[idx],color="anno",legend_loc="on data",save="_meta_endo_anno.png")

mk={"venules": ["AFF3","ACKR1","LRRC1","ADGRG6", "MMRN1", "POSTN"], "artery":["SEMA3G", "HEY1", "BMX"], "fenestrate": ["PLVAP", "CA4"], "postcap": ["PKHD1L1", "MMRN1", "ACKR1"], "arteirole": ["VEGFC", "NKAIN2"], "cap": ["SLC7A5", "BTNL9","INPP5D"] }

sc.pl.dotplot(adata_on, mk, groupby="leiden", save="_Endothelial_cell_hvg10000_epochnone_seurat_rs_0.6_clean_rm_sb_none_mk_all.png")

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_0.6_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad")

obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_0.5_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad.obs.gz", header=0, index_col=0)

adata.obs["leiden1"]=obs.loc[adata.obs.index,"leiden"].astype(str)


obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_1_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad.obs.gz", header=0, index_col=0)

idx=obs.loc[obs["leiden"].isin([3,9])].index
adata.obs.loc[idx,"leiden1"]="leiden"+"_"+obs.loc[idx,"leiden"].astype("str")

adata.obs.loc[adata.obs["leiden1"]=="3","leiden1"]="leiden_3"


sc.pl.umap(adata, color="leiden1", save="Endothelial_cell_hvg10000_epochnone_seurat_rs_0.5_1_clean_rm_sb_rm_4_6_none_scvi_trg.png")


sc.pl.dotplot(adata, mk, groupby="leiden1", save="_Endothelial_cell_hvg10000_epochnone_seurat_rs_0.5_1_clean_rm_sb_none_mk_all.png")


mk={"venules": ["VCAM1", "SELE", "ICAM1", "PECAM1", "CD34" ,"AFF3","ACKR1","LRRC1","ADGRG6", "MMRN1", "POSTN"], "artery":["SEMA3G", "HEY1", "BMX"], "fenestrate": ["PLVAP", "CA4"], "postcap": ["PKHD1L1", "MMRN1", "ACKR1"], "arteirole": ["VEGFC", "NKAIN2"], "cap": ["SLC7A5", "BTNL9","INPP5D"] }


mk={"venules": ["ACKR1","LRRC1","ADGRG6", "MMRN1", "POSTN"], "artery":["SEMA3G", "HEY1", "BMX"], "fenestrate": ["PLVAP", "CA4"], "postcap": ["PKHD1L1", "MMRN1", "ACKR1"], "arteirole": ["VEGFC", "NKAIN2"], "cap": ["SLC7A5", "BTNL9","INPP5D"] }


adata=sc.read("HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_v3_rs_1_subclass_sb_seurat_rmRPE_scvi_trg.h5ad")
bname="Endothelial_cell_hvg10000_epochnone_seurat_v3_rs_1.2_subclass_sb_seurat_rmRPE"

obs=pd.read_csv("HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_v3_rs_1.2_subclass_sb_seurat_rmRPE.obs.gz", header=0, index_col=0)
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/figures"
adata.obs["leiden1"]=obs.loc[adata.obs.index, "leiden"].astype(str)
sc.pl.dotplot(adata,mk,groupby="leiden1", save=f"_{bname}_mk1.png")

mk={"venules": ["NR2F2", "VCAM1", "ACKR1", "SELP",  "POSTN"], "artery":["EFNB2", "SOX17", "BMX", "SEMA3G", "HEY1",  "LTBP4", "FBLN5", "GJA5", "GJA4"], "fenestrate": ["PLVAP", "CA4"], "postcap": ["PKHD1L1", "MMRN1", "ACKR1"], "arteirole": ["VEGFC", "NKAIN2"], "cap": ["CA4", "PRX", "RGCC", "SPARC", "SGK1"] }

mk1={"artery":["CLDN10", "GJA5", "GJA4", "FBLIM1", "FBLN5", "FBLN2", "MGP", "BGN", "LTBP4", "LTBP1", "FN1", "SERPINE2", "CPAMD8", "CXCL12", "EFBN2", "SEMA3G", "VEGFA", "NOS1", "PDE3A", "PDE4D", "DKK2", "DLL4", "HEY1", "SOX5", "SOX17", "HES4", "PRDM16"]}
