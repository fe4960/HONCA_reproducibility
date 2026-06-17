import scanpy as sc

from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
rcParams.update({'font.size': 14})

#sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/"
#########adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/oligo_opc_seurat_new_scvi_hvg_2000_raw_velocity.h5ad")
##########desired_order = ["OLIGO2_LRRC7+", "OLIGO1", "OLIGO1_SVEP1+", "OLIGO1_RBFOX1+"]
########sc.pl.violin(adata[adata.obs["majorclass"]=="Oligodendrocyte"],"latent_time",groupby="subclass", order=desired_order, save="_opc_latent_time.pdf",rotation=90, jitter=False)


##########sc.pl.violin(adata[adata.obs["majorclass"]=="Oligodendrocyte"],"latent_time",groupby="subclass", order=desired_order, save="_opc_latent_time.png")
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/veloc/"

#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/veloc/Oligodendrocyte_precursor_cell_subclass_seurat_hvg_2000_raw_velocity_reservse.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/veloc/Oligodendrocyte_precursor_cell_subclass_seurat_cycling_hvg_2000_raw_velocity_reservse.h5ad")
adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]
dt={"COP_LAMA2+":"NFOL_UTRN+"}
adata.obs["subclass"]=adata.obs["subclass"].replace(dt)
order=["OPC_GLIS3+","OPC","OPC_Cycling", "COP_NFOL", "NFOL_UTRN+", "MFOL"]
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/veloc/"
sc.pl.violin(adata, "latent_time",groupby="subclass", order=order, save="_opc2_latent_time1.pdf",rotation=90, jitter=False)
