import scanpy as sc
from matplotlib import rcParams
import pandas as pd
rcParams["figure.figsize"] = (5,5)
sc.settings.figdir = "HCA_ON/data/5_refine_major/scvi/major/figures"
#fl=["hvg10000_epochnone_seurat_v3_rs","hvg2000_epoch20_seurat_v3_rs","hvg2000_epochnone_seurat_rs","hvg2000_epochnone_seurat_v3_rs"]
#fl=["hvg10000_epoch20_seurat_v3_rs"]
#fl=["hvg10000_epoch20_seurat_rs"]
#fl=["hvg2000_nonsb_epoch20_seurat_v3_rs","hvg10000_sb_epochnone_seurat_v3_rs"] #, "hvg2000_nonsb_epochnone_seurat_v3_rs"]
fl=["hvg10000_sb_epochnone_seurat_v3_rs"]
obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/ON_ONH_clean_new.obs.gz",index_col=0,header=0)
for f1 in fl:
    adata=sc.read(f"HCA_ON/data/5_refine_major/scvi/major/major_{f1}_1_clean_ON_ONH_new_scvi_trg.h5ad")
#    adata.obs["subclass"]=obs.loc[adata.obs.index,"subclass"]
    #    adata=sc.read(f"HCA_ON/data/5_refine_major/scvi/major/major_{f1}_1_clean_ON_ONH_scvi_trg.h5ad")
#    adata=sc.read(f"HCA_ON/data/5_refine_major/scvi/major/major_{f1}_1_clean_ON_ONH_sb_scvi_trg.h5ad")
   
#    adata.obs["majorclass1"]=adata.obs["majorclass1"].cat.add_categories(["Adipocyte", "Pigmented_cell", "Unknown?", "T_cell", "NK_cell","Mast_cell","Dendritic_cell","B_cell_activated"])
#    adata.obs["majorclass1"]=adata.obs["majorclass1"].cat.add_categories([ "Pigmented_cell", "Unknown?", "T_cell", "NK_cell","Mast_cell","Dendritic_cell","B_cell_activated","B_cell"])

#    adata.obs.loc[adata.obs["subclass"]=="Adipocyte","majorclass1"]="Adipocyte"
#    adata.obs.loc[adata.obs["subclass"].isin(["Pigmented_cell_1","Pigmented_cell_2"]),"majorclass1"]="Pigmented_cell"
    
#    adata.obs.loc[adata.obs["majorclass1"]=="NK_T_cell","majorclass1"]=adata.obs.loc[adata.obs["majorclass1"]=="NK_T_cell","subclass"].astype(str)

#    adata.obs.loc[adata.obs["majorclass1"]=="B_cell_activated","majorclass1"]="B_cell"

#    adata.obs["majorclass1"]=adata.obs["majorclass1"].cat.remove_unused_categories()

  # sc.pl.embedding(adata, basis="X_umap", color=["majorclass1"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_majorclass.png")
  # sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_subclass.png")

    sc.pl.embedding(adata, basis="X_umap", color=["majorclass1"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_new_majorclass.png")
#    sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_new_subclass.png")


#    sc.pl.embedding(adata, basis="X_umap", color=["majorclass1"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_sb_majorclass.png")
#    sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_sb_subclass.png")

