import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_raw_normcount_only.h5ad")
adata.obs["celltype"].value_counts()
adata=adata[adata.obs["celltype"]!="Unknown?"]
adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_raw_normcount_only_rmUnk.h5ad")
