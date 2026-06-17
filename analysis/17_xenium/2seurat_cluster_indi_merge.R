options(Seurat.object.assay.version = 'v3')
library(Seurat)
library(ggplot2)
library(dplyr)
#library(arrow)
args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH" #args[2]
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
#md=matching_dirs


xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_cutoff20.rds"))
label=args[2]

#####i=1

######bc=read.table(file=paste0(args[1],"_",i),header=T,row.names=1)
######xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)
#xenium_data_list[[i]]=xenium_data_list[[i]][features,]

######merged_obj=xenium_data_list[[i]]

for(i in 1:length(xenium_data_list)){

bc=read.table(file=paste0(args[1],"_",i),header=T,row.names=1)
xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)
#xenium_data_list[[i]]=xenium_data_list[[i]][features,]
#merged_obj <- merge(x = merged_obj, y = xenium_data_list[[i]] )


#new_assay <- CreateAssayObject(counts = xenium_data_list[[i]]@assays[["Xenium"]]$count)

#xenium_data_list[[i]]@assays[["RNA"]]=new_assay
#Key(xenium_data_list[[i]]@assays[["RNA"]])="RNA_"

#xenium_data_list[[i]]@assays[["Xenium"]]=NULL
####xenium_data_list[[i]]@assays[["BlankCodeword"]]=NULL
#####xenium_data_list[[i]]@assays[["ControlCodeword"]]=NULL
#####xenium_data_list[[i]]@assays[["ControlProbe"]]=NULL
#####xenium_data_list[[i]]@images[["fov"]]=NULL


seurat_obj <- CreateSeuratObject(counts = xenium_data_list[[i]]@assays[["Xenium"]]$counts)
seurat_obj@meta.data$sampleid=paste0("PP_",i)

#DefaultAssay(xenium_data_list[[i]]) <- "RNA"
#xenium_data_list[[i]]@meta.data$sampleid=paste0("PP_",i)
#saveRDS(xenium_data_list[[i]], paste0(output_PATH, target_tissue, "_xenium_merge_cutoff20_",label,"_",i,".rds"))
saveRDS(seurat_obj, paste0(output_PATH, target_tissue, "_xenium_merge_cutoff20_",label,"_",i,".rds"))

}


#saveRDS(merged_obj, paste0(output_PATH, target_tissue, "_xenium_merge_cutoff20_",label,".rds"))


