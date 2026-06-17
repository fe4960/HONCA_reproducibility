library(Seurat)
library(ggplot2)
library(dplyr)
#library(pheatmap)
#library(arrow)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_gap_stat_niches10_clean_merged_major.rds"))

#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_gap_stat_niches10_clean_merged.rds"))
id=c("PP_1", "PP_2", "ONH_1", "ONH_2", "ONH_3")

for(i in 1:length(id)){
x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]

fov="fov"
if(i>1){
fov=paste0("fov.",i)
}

DefaultBoundary(x1[[fov]]) <-  "segmentation"
#kmean7


pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche_clean_anno_k8_w1_major.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "kmean8", border.color=NA, axes=TRUE, border.size = 0, dark.background = T, cols = "Set1", fov=fov)
#cols = "polychrome"
print(p)
dev.off()



pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche_clean_anno_k6_w1_major.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "kmean6", border.color=NA, axes=TRUE, border.size = 0, dark.background = T, cols = "Set1", fov=fov)
#cols = "polychrome"
print(p)
dev.off()


pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche_clean_anno_k16_w1_major.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "kmean16", border.color=NA, axes=TRUE, border.size = 0, dark.background = T, cols = "polychrome", fov=fov)
#cols = "polychrome"
print(p)
dev.off()

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche_clean_anno_k10_w1_major.pdf"), width = 50, height = 35)

###pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche_clean_anno_harmony.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "kmean10", border.color=NA, axes=TRUE, border.size = 0, dark.background = T, cols = "polychrome", fov=fov)
#cols = "polychrome"
print(p)
dev.off()


pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche_clean_anno_k12_w1_major.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "kmean12", border.color=NA, axes=TRUE, border.size = 0, dark.background = T, cols = "polychrome", fov=fov)
print(p)
dev.off()

}
