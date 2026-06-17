library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(factoextra)
#library(arrow)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
############

cut=20

########
#####process xenium data

################
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")



seu=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_gap_stat_niches10_clean_merged_major.rds"))
library(viridis)

dt=table(seu$harmony_anno, seu$kmean16)

colnames(dt)=seq(1,16,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "row",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = T,    # disable column clustering
  main = paste0("kmean=16"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean16_clean_major_scale_row.pdf"), width=5, height=4)
print(p)
dev.off()


dt=table(seu$harmony_anno, seu$kmean6)

colnames(dt)=seq(1,6,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "row",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=6"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean6_clean_major_scale_row.pdf"), width=5, height=3)
print(p)
dev.off()

dt=table(seu$harmony_anno, seu$kmean10)

colnames(dt)=seq(1,10,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "row",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=10"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean10_clean_major_scale_row.pdf"), width=5, height=3)
print(p)
dev.off()

dt=table(seu$harmony_anno, seu$kmean12)

colnames(dt)=seq(1,12,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "row",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=12"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean12_clean_major_scale_row.pdf"), width=5, height=3)
print(p)
dev.off()


dt=table(seu$harmony_anno, seu$kmean8)

colnames(dt)=seq(1,8,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "row",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=8"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean8_clean_major_scale_row.pdf"), width=5, height=3)
print(p)
dev.off()

