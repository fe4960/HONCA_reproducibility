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

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))

obs=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/sn_st_scvi_harmony_obs.txt.gz", sep="\t", header=T, row.names=1)

score=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/sn_st_scvi_harmony_score.txt.gz", sep="\t", header=T, row.names=1)
row_max <- apply(score, 1, max)

obs=obs[row_max >=0.7,]

library(ggplot2)

for(i in 1:length(xenium_data_list)){

xenium_obj=xenium_data_list[[i]]

slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

#merged_obj1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

obs1=obs[obs$sampleid==slide_id,]
#meta=x1@meta.data #merged_obj@meta.data[merged_obj@meta.data$slide_id==id[i],]

xenium_obj=xenium_obj[,obs1$barcode]

xenium_obj@meta.data$lr_celltype="other"
#xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$harmony_anno=merged_obj1@meta.data$harmony_anno
xenium_obj@meta.data[obs1$barcode,]$lr_celltype=obs1$lr_celltype
xenium_obj=xenium_obj[,!(xenium_obj@meta.data$lr_celltype %in% c("Pigmented_cell","Adipocyte"))]

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", slide_id, "_scvi_harmony.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "lr_celltype", border.size = 0,   cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )

#+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(p)
dev.off()

}

#####

#========

