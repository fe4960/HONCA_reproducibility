import scanpy as sc

#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_normcount.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount.h5ad")

cell=["Astrocyte", "Endothelial_cell", "Fibroblast", "Microglia", "Mural_cell", "Oligodendrocyte", "Oligodendrocyte_precursor_cell"]

adata=adata[adata.obs["majorclass"].isin(cell)]


adata.obs["subclass"]=adata.obs["majorclass"].astype(str)

idx=adata.obs.loc[adata.obs["majorclass"].isin(["Astrocyte","Oligodendrocyte"])].index

adata.obs.loc[idx,"subclass"]=adata.obs.loc[idx,"celltype"].astype(str)


adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.h5ad")
adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.obs.gz")

#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_normcount_subset.h5ad")

#adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_normcount_subset.obs.gz")


######
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount.h5ad")

cell=["Astrocyte", "Endothelial_cell", "Fibroblast", "Microglia", "Mural_cell", "Oligodendrocyte", "Oligodendrocyte_precursor_cell"]

adata=adata[adata.obs["majorclass"].isin(cell)]


adata.obs["subclass"]=adata.obs["majorclass"].astype(str)

idx=adata.obs.loc[adata.obs["majorclass"].isin(["Astrocyte","Oligodendrocyte", "Fibroblast"])].index

adata.obs.loc[idx,"subclass"]=adata.obs.loc[idx,"celltype"].astype(str)


adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.h5ad")
adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.obs.gz")
