adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount_rename_majorclass.h5ad")
adata.obs["majorclass1"]=adata.obs["majorclass1"].str.replace("/","_")
adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount_rename_majorclass_final.h5ad")
