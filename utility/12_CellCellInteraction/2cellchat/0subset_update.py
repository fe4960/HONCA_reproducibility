import scanpy as sc
import pandas as pd


#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_raw_normcount_only.h5ad")

#cell=["Astrocyte", "Endothelial_cell", "Fibroblast", "Microglia", "Mural_cell", "Oligodendrocyte", "Oligodendrocyte_precursor_cell", "RPE"]
cell=["Astrocyte", "Endothelial_cell", "Fibroblast", "Microglia", "Mural_cell", "Oligodendrocyte", "Oligodendrocyte_precursor_cell"]

adata=adata[adata.obs["majorclass1"].isin(cell)]


adata.obs["subclass"]=adata.obs["majorclass1"].astype(str)

idx=adata.obs.loc[adata.obs["majorclass1"].isin(["Astrocyte","Oligodendrocyte", "Fibroblast"])].index

adata.obs.loc[idx,"subclass"]=adata.obs.loc[idx,"celltype"].astype(str)

adata.obs.loc[adata.obs["subclass"].str.contains("Astro_ON") & -adata.obs["subclass"].str.contains("Astro_ONH"), "subclass"]="Astro_ON"


adata.obs.loc[adata.obs["subclass"].str.contains("Astro_ONH") & -adata.obs["subclass"].str.contains("Astro_ONHON"), "subclass"]="Astro_ONH"

adata.obs.loc[adata.obs["subclass"].str.contains("Astro_retina") , "subclass"]="Astro_retina"


adata.obs.loc[adata.obs["subclass"].str.contains("OLIGO1") , "subclass"]="OLIGO1"

adata.obs.loc[adata.obs["subclass"].str.contains("arachnoid") , "subclass"]="Fibro_arachnoid"

adata.obs.loc[adata.obs["subclass"].isin(["Fibro_intima_pia","Fibro_epipial"]) , "subclass"]="Fibro_pia"

adata.obs.loc[adata.obs["subclass"].str.contains("Fibro_dura") , "subclass"]="Fibro_dura"

adata.X=adata.layers["counts"].copy()
del adata.layers
del adata.raw
#adata.obs.loc[adata.obs["subclass"].str.contains("RPE") , "subclass"]="RPE"

#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.h5ad")
#adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.obs.gz")
adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset.h5ad")
adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset.obs.gz")
