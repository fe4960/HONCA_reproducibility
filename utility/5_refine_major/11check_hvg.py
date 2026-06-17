import scanpy as sc
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg10000_epochnone_seurat_rs_1_clean_ON_ONH_scvi_trg.h5ad")
#adata.var["highly_variable"].to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg10000_epochnone_seurat_rs_1_clean_ON_ONH_scvi_trg_hvg.csv")

#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg10000_epochnone_seurat_v3_rs_1_clean_ON_ONH_scvi_trg.h5ad")
#adata.var["highly_variable"].to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg10000_epochnone_seurat_v3_rs_1_clean_ON_ONH_scvi_trg_hvg.csv")


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg2000_epochnone_seurat_rs_1_clean_ON_ONH_scvi_trg.h5ad")
adata.var["highly_variable"].to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg2000_epochnone_seurat_rs_1_clean_ON_ONH_scvi_trg_hvg.csv")

