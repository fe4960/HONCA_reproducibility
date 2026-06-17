import scanpy as sc
#sampleid=["BCM_22_0784_ON_RNA","MMD_23_17738_ON_RNA","BCM_23_0491_ONH_RNA","MMD_19_D008_ONH_RNA","MMD_23_22486_ONH_RNA","BCM_22_0890_ON_RNA"]
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_rawcount.h5ad")
#adata=adata[adata.obs["sampleid"].isin(sampleid)]
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad")
adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz")
