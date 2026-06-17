import scanpy as sc
import pandas as pd


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad")

mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_astro_subclass3_n_markers_rctd_top4.csv",sep=",",header=0, index_col=0)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"

adata.obs["subclass3"]="astro_retina"

adata.obs.loc[adata.obs["subclass"].str.contains("ONHON"), "subclass3"]="astro_ONHON"
adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ONH1"), "subclass3"]="astro_ONH"
adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ONH2"), "subclass3"]="astro_ONH"
adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ON") & ~adata.obs["subclass"].str.contains("ONH")  , "subclass3"]="astro_ON"

sc.pl.dotplot(adata, mk["gene"], groupby="subclass3", save= "_ONONH_comb_astro_subclass3_markers_rctd_ref.png" )
sc.pl.dotplot(adata, mk["gene"], groupby="subclass3", save= "_ONONH_comb_astro_subclass3_markers_rctd_ref_st_var.png", standard_scale="var" )


