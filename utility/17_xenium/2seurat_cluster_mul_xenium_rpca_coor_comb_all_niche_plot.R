library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)

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

cut=20

########
#####process xenium data

label="oligo"

################
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6", ".rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

id=na.omit(id)



for(i in 1:length(id)){
#x1= readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_", id[i], "_niches10.rds"))
#dt=read.table(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_", id[i], "_niches10.txt"), sep="\t" , header=T)

#########dt=read.table( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches10_clean_major.txt"), sep="\t" , header=T)

##########x1=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches10_clean_major.rds"))

dt=read.table( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_", id[i], "_niches10_clean_subclass.txt"),  sep="\t" , header=T)

#dt=read.table( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches10_clean.txt"), sep="\t" , header=T)

x1=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches10_clean.rds"))



fov="fov"
if(i>1){
fov=paste0("fov.",i)
}

dt=data.matrix(dt)

library(RColorBrewer)

paired_more = brewer.pal(12, "Paired")

nic <- as.character(seq(1, 10, 1))
#cols_named <- setNames(paired_more[1:10], nic)
cols_named <- setNames(paired_more[10:1], nic)


DefaultBoundary(x1[[fov]]) <-  "segmentation"

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche10_w1_sub.pdf"), width = 50, height = 35)

#########pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche10_w1.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "niches", size = 1.5, dark.background = T, cols = paired_more, fov=fov)
#p=ImageDimPlot(object = x1, group.by = "harmony_anno", border.size = 0,   cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )
#+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(p)
dev.off()
#n=ncol(dt)*nrow(dt)
colnames(dt)=seq(1,10,1)
n=100
#p=heatmap(
#  dt,
#  col = viridis(n),    # viridis color palette
#  scale = "none",        # don't rescale rows or columns
#  margins = c(5, 5),     # adjust margin size
#  main = paste0("majorclass_",id[i]),
#  cexRow = 1,
#  cexCol = 1
#)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("major_", id[i]),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)


pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_heatmap_niche10_w1_sub.pdf"), width=5, height=5)

###########pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_heatmap_niche10_w1.pdf"), width=5, height=5)
print(p)
dev.off()

dt=dt[!(rownames(dt) %in% c("Rod","Cone","BC","HC","MG","RGC","AC", "Melanocyte","Schwann_cell")), ]

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("major_", id[i]),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_heatmap_foc_niche10_w1_sub.pdf"), width=5, height=3)

###########pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_heatmap_foc_niche10_w1.pdf"), width=5, height=3)
print(p)
dev.off()




}
