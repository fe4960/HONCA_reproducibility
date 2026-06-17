import scanpy as sc
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)


adata=sc.read("HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location1.h5ad")
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"

mk=["ACTA2", "TPM1"]
sc.pl.dotplot(adata, mk, groupby="subclass", save=f"_subclass_Fibro_HNMN1.pdf", swap_axes=True)


