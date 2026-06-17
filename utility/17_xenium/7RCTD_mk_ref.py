import scanpy as sc
import pandas as pd

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat.h5ad")
#fn=
#mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_oligo_subclass_015_markers_rctd.csv",sep=",",header=0, index_col=0)
mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_oligo_subclass_01_markers_rctd.csv",sep=",",header=0, index_col=0)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"

#sc.pl.dotplot(adata, mk["gene"], groupby="subclass", save= "_ONONH_comb_oligo_subclass_015_markers_rctd_ref.png" )
sc.pl.dotplot(adata, mk["gene"], groupby="subclass", save= "_ONONH_comb_oligo_subclass_01_markers_rctd_ref.png" )

adata.obs["subclass1"]=adata.obs["subclass"]

adata.obs.loc[(adata.obs["subclass"]=="OLIGO1_RBFOX1+") | (adata.obs["subclass"]=="OLIGO1_SVEP1+"), "subclass1"] = "OLIGO1"


adata.obs["subclass1"] = adata.obs["subclass1"].cat.remove_unused_categories()

#sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_oligo_subclass_015_markers_rctd_ref_two_class.png" )
sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_oligo_subclass_01_markers_rctd_ref_two_class.png" )

#sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_oligo_subclass_015_markers_rctd_ref_two_class_st_var.png", standard_scale="var" )
sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_oligo_subclass_01_markers_rctd_ref_two_class_st_var.png", standard_scale="var" )


mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_oligo_subclass1__01_markers_rctd.csv",sep=",",header=0, index_col=0)

sc.pl.dotplot(adata, mk["gene"], groupby="subclass", save= "_ONONH_comb_oligo_subclass1_01_markers_rctd_ref.png" )
sc.pl.dotplot(adata, mk["gene"], groupby="subclass", save= "_ONONH_comb_oligo_subclass1_01_markers_rctd_ref_st_var.png", standard_scale="var" )

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad")

#mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_astro_subclass1_n_markers_rctd_top5.csv",sep=",",header=0, index_col=0)
mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_astro_subclass1_n_markers_rctd_top4.csv",sep=",",header=0, index_col=0)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"

adata.obs["subclass1"]="astro_retina"

adata.obs.loc[adata.obs["subclass"].str.contains("ONHON"), "subclass1"]="astro_ONHON"
adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ONH1"), "subclass1"]="astro_ONH1"
adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ONH2"), "subclass1"]="astro_ONH2"
adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ON") & ~adata.obs["subclass"].str.contains("ONH")  , "subclass1"]="astro_ON"


#data$subclass1[grepl("retina",data$subclass)]="astro_retina"
#data$subclass1[grepl("ONHON",data$subclass)]="astro_ONHON"
#data$subclass1[(!grepl("ONHON",data$subclass))& grepl("ONH1",data$subclass)]="astro_ONH1"
#data$subclass1[(!grepl("ONHON",data$subclass))& grepl("ONH2",data$subclass)]="astro_ONH2"
#data$subclass1[(!grepl("ONHON",data$subclass))& grepl("ON",data$subclass)& ( !grepl("ONH",data$subclass))]="astro_ON"



sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_astro_subclass1_markers_rctd_ref.png" )

sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_astro_subclass1_markers_rctd_ref_st_var.png", standard_scale="var" )


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new.h5ad")

mk=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_fibro_subclass1_n_markers_rctd_top4.csv",sep=",",header=0, index_col=0)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"

adata.obs["subclass1"]=[ x.split("_")[0]+"_"+ x.split("_")[1] for x in adata.obs["subclass"]]

#adata.obs.loc[adata.obs["subclass"].str.contains("ONHON"), "subclass1"]="astro_ONHON"
#adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ONH1"), "subclass1"]="astro_ONH1"
#adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ONH2"), "subclass1"]="astro_ONH2"
#adata.obs.loc[~adata.obs["subclass"].str.contains("ONHON") & adata.obs["subclass"].str.contains("ON") & ~adata.obs["subclass"].str.contains("ONH")  , "subclass1"]="astro_ON"

sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_fibro_subclass1_markers_rctd_ref.png" )

sc.pl.dotplot(adata, mk["gene"], groupby="subclass1", save= "_ONONH_comb_fibro_subclass1_markers_rctd_ref_st_var.png", standard_scale="var" )
