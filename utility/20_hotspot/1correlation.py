import pandas as pd
import scanpy as sc
import numpy as np
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/Oligodendrocyte_subclass_seurat_hvg_2000_raw_velocity.h5ad"


adata=sc.read(h5ad)

cor=adata.obs["age"].corr(adata.obs["velocity_pseudotime"])
#0.05270479727365323

cor=adata.obs["age"].corr(adata.obs["latent_time"])
#0.034712280793621084

