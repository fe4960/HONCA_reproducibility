import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
celltype="Oligodendrocyte"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
mk="sanes_mk"
#adata=sc.read_h5ad(f"{outdir}/{celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
#HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_hvg10000_epochnone_seurat_rs_0.4_clean_sb_none_scvi_trg.h5ad
adata1=sc.read_h5ad(f"{outdir}/oligo_opc_hvg10000_epochnone_seurat_rs_0.4_clean_sb_none_scvi_trg.h5ad")
adata=sc.read_h5ad(f"{outdir}/oligo_opc_seurat_new.h5ad")

obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_hvg10000_epochnone_seurat_rs_0.8_clean_sb_none.obs.gz",header=0, index_col=0)
idx=obs.loc[obs["leiden"]==10].index


adata.obs["subclass"]=adata.obs["subclass"].cat.add_categories(["OPC_Cycling"])
adata.obs.loc[idx,"subclass"]="OPC_Cycling"

bname=f"oligo_opc_seurat_new_scvi_cycling"

#####bname=f"oligo_opc_seurat_new_scvi"
#adata.obs["subclass"]="Oligodendrocyte_precursor_cell"
sc.settings.figdir = f"{outdir}/figures"

dt={"OLIGO_LRRC7+":"OLIGO2_LRRC7+", "OLIGO_non_specific":"OLIGO1", "OLIGO_SVEP1+":"OLIGO1_SVEP1+", "OLIGO_RBFOX1+":"OLIGO1_RBFOX1+"}
adata.obs["subclass"]=adata.obs["subclass"].replace(dt)
adata=adata[adata1.obs.index]




del adata.obsm["X_scVI"]
adata.obsm["X_scVI"]=adata1.obsm["X_scVI"]
del adata.obsm["X_umap"]
adata.obsm["X_umap"]=adata1.obsm["X_umap"]
##del adata.uns["umap"]
adata.uns["umap"]=adata1.uns["umap"]
del adata.obs["leiden"]
adata.obs["leiden"]=adata1.obs["leiden"]

adata.write(f"{outdir}/{bname}.h5ad")
sc.pl.umap(adata,color="subclass",save="_opc_oligo_scvi.png", palette="tab20")
