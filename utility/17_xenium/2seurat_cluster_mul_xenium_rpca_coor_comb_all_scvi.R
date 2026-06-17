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


xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))

#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6", ".rds"))

#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

#id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

#id=na.omit(id)

for(i in 1:5){
#x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]
x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]

#x1@meta.data$x_coor=x1@images$fov$centroids$x


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
#x1_meta=x1_meta[c("AC","Astrocyte","BC","Cone","Endothelial-cell","Fibroblast","HC","Immune-cell","Macrophage","MG","Microglia","Mural-cell","Oligodendrocyte","OPC","RGC","RPE"),]
#x1_meta=x1_meta[x1_meta$harmony_anno%in% c("AC","Astrocyte","BC","Cone","Endothelial_cell","Fibroblast","HC","Immune_cell","Macrophage","Melanocyte","MG","Microglia","Mural_cell","Oligodendrocyte","OPC","RGC","RPE", "Rod", "Schwann_cell"),]
#######x1_meta=x1_meta[x1_meta$harmony_anno%in% c("Astrocyte","BC","Cone","Endothelial_cell","Fibroblast","HC","Immune_cell","Macrophage","Melanocyte","MG","Microglia","Mural_cell","Oligodendrocyte","OPC","RGC","RPE", "Rod", "Schwann_cell"),]

x1_meta_all=rbind(x1_meta_all,x1_meta)
}


#x1_meta_all[x1_meta_all$xenium_clusters_p2==4,]$xenium_clusters_p2=3
#x1_meta_all[x1_meta_all$xenium_clusters_p2==5,]$xenium_clusters_p2=4
#x1_meta_all[x1_meta_all$xenium_clusters_p2==6,]$xenium_clusters_p2=0


#x1_meta_all=x1_meta_all[c("AC","Astrocyte","BC","Cone","Endothelial-cell","Fibroblast","HC","Immune-cell","Macrophage","MG","Microglia","Mural-cell","Oligodendrocyte","OPC","RGC","RPE"),]


df <- x1_meta_all %>%
  group_by(bin_y, harmony_anno) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


#ggplot(df,aes(x=bin,y=prop,fill=xenium_clusters_p2))+geom_bar(stat = "identity") +
#  ylab("Proportion") +
#  theme_minimal()
library(RColorBrewer)

#paired_more <- colorRampPalette(brewer.pal(8, "Set2"))(19)

paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))
#paired_more=paired_more[-1]
#paired_more=paired_more[-1]
#paired_more=paired_more[3:]
fn=paste0(output_PATH, target_tissue,"_cell_cluster_dist_y_comb_p6.pdf")
p=ggplot(df,aes(x=bin_y,y=prop,fill=harmony_anno))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_manual(values=paired_more)
pdf(fn)
print(p)
dev.off()

#########


#for(i in 1:length(id)){

#}




#}
#summary(lm(df[df$xenium_clusters_p2==0,]$prop~df[df$xenium_clusters_p2==0,]$bin))
