import scanpy as sc
import anndata as ad
mic=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/Microglia_subclass_sb.h5ad")
#mic.X=mic.layers["counts"]

mac=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/Macrophage_subclass_sb.h5ad")
#mac.X=mac.layers["counts"]



adata_list=[mic, mac]
adata=ad.concat(adata_list)

adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/Micro_Macro_clean_final.h5ad")
