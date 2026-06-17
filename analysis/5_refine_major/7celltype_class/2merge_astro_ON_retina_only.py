import scanpy as sc
import anndata as ad
import pandas as pd
#adata=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HRCA/Lattice/Retina_h5ad/HRCA_snRNA_allcells_rawcounts.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HRCA/CAP/fix_donor/scrnah5adobssubset2replacebydict/HRCA_snRNA_allcells_rawcounts.h5ad")
ret_astro=adata[adata.obs["majorclass"]=="Astrocyte"]
ret_astro.obs["subclass"]="Astro_retina_HRCA"


#ret_mg=adata[adata.obs["majorclass"]=="MG"]
#ret_mg.obs["subclass"]="MG_retina_HRCA"


#on_astro=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new.h5ad")
on_astro=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean.h5ad")

on_astro.X=on_astro.layers["counts"]


#mg=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/MG/clean/MG_subclass_sb.h5ad")
#mg=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/MG/clean/MG_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
#obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/MG/clean/MG_hvg10000_epochnone_seurat_rs_0.4_rmnone_sb.obs.gz",header=0,index_col=0)
#mg.X=mg.layers["counts"]
#mg.obs["subclass"]="MG_ONOH_"+obs.loc[mg.obs.index,"leiden"].astype("str")

astro=ad.concat([on_astro, ret_astro])

#astro.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astro_ret_ONONH.h5ad")

astro.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astro_ret_ONONH_update.h5ad")
