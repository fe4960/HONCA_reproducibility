# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)


od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="5" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue1, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)

cut=20
xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", ".rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

#id=unique(xenium.combined@meta.data$slide_id)
library(dplyr)
library(ggplot2)
#for(i in 2:length(id)){
id="xenium_ONH_1"
merged_obj1=xenium.combined[,xenium.combined@meta.data$slide==id]

cb=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/Selection_1_cells_stats_onh.csv",sep=",",header=T)

#merged_obj2=merged_obj1[, merged_obj1@meta.data$barcode%in%cb$Cell.ID]

rownames(cb)=cb[,"Cell.ID"]
merged_obj2 <- subset(merged_obj1, subset = barcode  %in% cb$Cell.ID)

cell=Cells(merged_obj1)
#Cells(merged_obj1)=merged_obj1@meta.data$barcode
#merged_obj2=merged_obj1[,cb$Cell.ID]

cb$Cell.ID=paste0(cb$Cell.ID,"_",3)

#merged_obj2 <- subset(merged_obj1, cells = cb$Cell.ID)

#merged_obj2=merged_obj1[,cb$Cell.ID]

df=data.frame(merged_obj1@meta.data)
df=df[,c("barcode","cell_size","slide","harmony_anno","slide_id")]
n=i-1
write.table(df,file=paste0(output_PATH, target_tissue, "coembed.combined_inte_", n, ".txt"),sep="\t",quote=F)

#}

