import scanpy as sc
dir="/dfs3b/ruic20_lab/junw42/"
adata=sc.read(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_normCount.h5ad")
cell=["Muller glia"]
for c in cell:
    name=c.replace(" ","_")
    adata[adata.obs["majorclass"]==c].write(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/{name}.h5ad")
