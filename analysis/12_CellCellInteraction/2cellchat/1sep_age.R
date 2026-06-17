library(Seurat)
library(dplyr)
#lattice_update_final/ONONH_all_rawcount_subset_simple
#data=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_subset.rds")
data=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset.rds")

data@meta.data$subclass1=data@meta.data$subclass

data@meta.data$subclass1[grepl("Astro_ON", data@meta.data$subclass) & (!grepl("Astro_ONHON", data@meta.data$subclass)) & (!grepl("Astro_ONH", data@meta.data$subclass))] = "Astro_ON"

data@meta.data$subclass1[ (grepl("ONHON", data@meta.data$subclass))] = "Astro_ONHON"

data@meta.data$subclass1[grepl("retina", data@meta.data$subclass)] = "Astro_retina"

data@meta.data$subclass1[grepl("OLIGO1", data@meta.data$subclass)] = "OLIGO1"

data@meta.data$subclass1[data@meta.data$subclass=="Fibroblast"]=data@meta.data$celltype[data@meta.data$subclass=="Fibroblast"]


data@meta.data$subclass1[grepl("Fibro_pia", data@meta.data$celltype)]="Fibro_pia"

data@meta.data$subclass1[grepl("Fibro_sclera", data@meta.data$celltype)]="Fibro_sclera"  

data@meta.data$subclass1[grepl("Fibro_arachnoid", data@meta.data$celltype)]="Fibro_arachnoid"  



data1=subset(data, subset = age_year <= 12)

saveRDS(data1, "/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset_age0to12.rds")



data1=subset(data, subset = age_year >= 20 & age_year <= 40 )


saveRDS(data1, "/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset_age20to40.rds")


data1=subset(data, subset = age_year >= 70 )



saveRDS(data1, "/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset_age70.rds")
