import scanpy as sc
data=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount.h5ad")
data.obs["majorclass1"]=data.obs["majorclass"].str.replace(" ", "_")
del data.uns
del data.obsm
del data.obsp
data.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount_rename_majorclass.h5ad")
