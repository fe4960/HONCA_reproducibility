import scanpy as sc
import anndata as ad
#adata=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HRCA/Lattice/Retina_h5ad/HRCA_snRNA_allcells_rawcounts.h5ad")

#ret_astro=adata[adata.obs["majorclass"]=="Astrocyte"]
#ret_astro.obs["subclass"]="Astro_retina_HRCA"


#ret_mg=adata[adata.obs["majorclass"]=="MG"]
#ret_mg.obs["subclass"]="MG_retina_HRCA"


on_astro=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new.h5ad")
on_astro=on_astro[on_astro.obs["subclass"].isin(["Astro_ON1","Astro_ON2","Astro_ON3"])]
#on_astro.X=on_astro.layers["counts"]


#mg=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/MG/clean/MG_subclass.h5ad")

#mg.X=mg.layers["counts"]



#astro=ad.concat([on_astro, ret_astro, mg, ret_mg])

#astro.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astro_ret_ONONH.h5ad")

on_astro.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte_ON/Astro_ON.h5ad")
