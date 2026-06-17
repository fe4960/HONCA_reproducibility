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
# read xenium data
#path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0)

#xenium.combined <- args[1]
#parts <- args[2]

#============================================================================================================
# co-embedding - Harmony


library(harmony)
library(RColorBrewer)




#merged_obj=readRDS(paste0(output_PATH, target_tissue, "_6_coembed.combined_cutoff20_ref_xenium_no_inte_harmony.rds"))
	


xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_6_coembed.combined_inte_harmony_w1", ".rds"))

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_6_xenium_data_list_comb_cutoff20.rds"))

cut=20

#xenium.combined=merged_obj   #readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1", ".rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

#id=unique(xenium.combined@meta.data$slide_id)

library(ggplot2)
#for(i in 1:(length(id)-1)){
#for(i in 1:(length(id)){
for(i in 1:length(xenium_data_list)){

xenium_obj=xenium_data_list[[i]]

slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

merged_obj1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]


#meta=x1@meta.data #merged_obj@meta.data[merged_obj@meta.data$slide_id==id[i],]

xenium_obj@meta.data$harmony_anno="other"
xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$harmony_anno=merged_obj1@meta.data$harmony_anno

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

paired_more = c("azure4", brewer.pal(12, "Set3"), "darkviolet", brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))

#paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))

celltypes <- sort(unique(xenium_obj$harmony_anno))

cols_named <- setNames(paired_more[1:length(celltypes)], celltypes)

pdf(paste0(output_PATH, target_tissue, "_6_all_annot_majorclass_", slide_id, "_w1_1.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0,   cols = cols_named, border.color=NA, axes=TRUE, dark.background=TRUE )
#"polychrome"
#+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(p)
dev.off()

}

#####

#========

