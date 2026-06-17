library(Seurat)
library(ggplot2)
library(dplyr)
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
cut=50



#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1", ".rds"))

xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_gap_stat_niches10_clean_merged_major.rds"))
#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_gap_stat_niches10_clean_merged.rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

id=na.omit(id)

for(i in 1:length(id)){
x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]


x1@meta.data$y_coor=x1@images$fov$centroids$y


if((i==2) || (i==5)){
x1@meta.data$y_coor=-x1@images$fov$centroids$y
}

#library(ggplot2)
############



x1@meta.data$bin_y <- cut(
  x1@meta.data$y_coor,
  breaks = 100,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE
)

x1_meta=x1@meta.data
x1_meta_all=rbind(x1_meta_all,x1_meta)
}




df <- x1_meta_all %>%
  group_by(bin_y, kmean12) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


library(RColorBrewer)


paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))
########fn=paste0(output_PATH, target_tissue,"_cell_cluster_dist_y_comb_w1_kmean12.pdf")
fn=paste0(output_PATH, target_tissue,"_cell_cluster_dist_y_comb_w1_kmean12_sub.pdf")

p=ggplot(df,aes(x=bin_y,y=prop,fill=as.factor(kmean12)))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_manual(values=paired_more)
pdf(fn)
print(p)
dev.off()


df <- x1_meta_all %>%
  group_by(bin_y, kmean10) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


library(RColorBrewer)


paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))
#########fn=paste0(output_PATH, target_tissue,"_cell_cluster_dist_y_comb_w1_kmean10.pdf")
fn=paste0(output_PATH, target_tissue,"_cell_cluster_dist_y_comb_w1_kmean10_sub.pdf")

paired_more = brewer.pal(10, "Paired")

p=ggplot(df,aes(x=bin_y,y=prop,fill=factor(kmean10)))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_manual(values=paired_more[10:1])
  #scale_fill_brewer(palette="Paired")
	#  ylab("Proportion") + theme_minimal() + scale_fill_manual(values=paired_more)
pdf(fn, height=6, width=6)
print(p)
dev.off()
