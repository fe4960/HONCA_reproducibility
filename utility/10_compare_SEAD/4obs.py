import scanpy as sc
import pandas as pd
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/"

cell=["Astrocytes", "Microglia", "OPCs"]

for c in cell:
    h5ad=f"{dir1}/GSE267301_{c}_Subset.h5ad"
    adata=sc.read(h5ad)
    adata.obs.to_csv(f"{dir1}/GSE267301_{c}_Subset.obs.gz")
