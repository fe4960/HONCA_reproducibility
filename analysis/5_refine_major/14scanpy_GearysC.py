import scanpy as sc
import sys
import numpy as np

main="/dfs3b/ruic20_lab/junw42"
#cell=sys.argv[1]
#h5ad=f"{main}/HCA_ON/data/5_refine_major/scvi/{cell}/clean/{cell}_hvg10000_epochnone_seurat_v3_rs_1_span_1_subclass_scvi_trg.h5ad"
#h5ad=f"{main}//HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_sb_scvi_trg.h5ad"
##########h5ad=f"{main}//HCA_ON/data/5_refine_major/scvi/major/major_hvg2000_nonsb_epoch20_seurat_v3_rs_1_clean_ON_ONH_new_scvi_trg.h5ad"
##########data=sc.read(h5ad)
#########print("hvg2000_nonsb_epoch20_seurat_v3_rs_1")
##########pc_c = sc.metrics.gearys_c(data, obsm="X_scVI")
##########print(pc_c)
#h5ad=f"{main}/HCA_ON/data/5_refine_major/scvi/{cell}/clean/{cell}_subclass.h5ad"
#h5ad=f"{main}//HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_scvi_trg.h5ad"
h5ad=f"{main}//HCA_ON/data/5_refine_major/scvi/major/major_hvg10000_sb_epochnone_seurat_v3_rs_1_clean_ON_ONH_new_scvi_trg.h5ad"
data=sc.read(h5ad)
print("hvg10000_sb_epochnone_seurat_v3_rs_1")
pc_c = sc.metrics.gearys_c(data, obsm="X_scVI")
print(pc_c)

h5ad=f"{main}//HCA_ON/data/5_refine_major/scvi/major/major_hvg2000_nonsb_epochnone_seurat_v3_rs_1_clean_ON_ONH_new_scvi_trg.h5ad"
data=sc.read(h5ad)
print("hvg2000_nonsb_epochnone_seurat_v3_rs_1")
pc_c = sc.metrics.gearys_c(data, obsm="X_scVI")
print(pc_c)
