import matplotlib.pyplot as plt
import scanpy as sc

#sys.setrecursionlimit(10000)

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

adata=sc.read("HCA_ON/data/5_refine_major/scvi/major/clean/major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_ret_scvi_trg.h5ad")

sc.settings.figdir="HCA_ON/data/5_refine_major/scvi/major/clean/figures/"

sc.pl.umap(adata, color="majorclass1", save="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_ret_majorclass.png", frameon=False)


sc.pl.umap(adata, color="majorclass1", save="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_ret_onData_majorclass.png", frameon=False, legend_loc="on data" )

