options(Seurat.object.assay.version = 'v3')
library(Seurat)
library(ggplot2)
library(dplyr)
#library(arrow)
args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="3" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
#target_tissue="ONONH"
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0) 
matching_dirs0 <- list.dirs(path, full.names = TRUE, recursive = FALSE)
xenium_data_list <- list()
matching_dirs=matching_dirs0[grepl(kw,matching_dirs0)]


xenium_data=readRDS(paste0(output_PATH, target_tissue,"_fibro_xenium.combined_cutoff20_fc_5.rds"))

#HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_fibro_xenium.combined_cutoff20_fc_5.rds
#xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))
#label=args[2]

id=unique(xenium_data@meta.data$slide_id)

#for(i in 1:length(xenium_data_list)){
for(i in id){
#bc=read.table(file=paste0(args[1],"_",i),header=T,row.names=1)
#xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)


#seurat_obj <- CreateSeuratObject(counts = xenium_data_list[[i]]@assays[["Xenium"]]$counts)

xenium_data1=subset(xenium_data, subset = slide_id == i)
seurat_obj <- CreateSeuratObject(counts = xenium_data1@assays[["Xenium"]]$counts)

seurat_obj@meta.data$sampleid=i #names(xenium_data_list)[i]

seurat_obj@meta.data$barcode=xenium_data1@meta.data$barcode 

saveRDS(seurat_obj, paste0(output_PATH, target_tissue, "_fibro_xenium.combined_cutoff20_fc_5_",i,".rds"))


#sam=paste0(target_tissue, "_fibro_xenium.combined_cutoff20_fc_5_",i)

#obj=readRDS(paste0(output_PATH, target_tissue, "_fibro_xenium.combined_cutoff20_fc_5_",i,".rds"))
#filename=paste0(output_PATH, target_tissue, "_fibro_xenium.combined_cutoff20_fc_5_",i,".h5Seurat") #output_PATH, target_tissue, "_coembed.combined_inte.h5Seurat")



}


