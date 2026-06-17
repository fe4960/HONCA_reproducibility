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

cut=50

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


xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)

#meta=xenium.combined@meta.data

#df=meta %>% group_by(xenium_clusters_p2) %>% summarize (mean=mean(nFeature_Xenium))

for(i in 1:length(id)){
x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]

x1@meta.data$x_coor=x1@images$fov$centroids$x
x1@meta.data$y_coor=x1@images$fov$centroids$y


#library(ggplot2)
x1@meta.data$bin_x <- cut(
  x1@meta.data$x_coor,
  breaks = 100,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE
)

x1_meta=x1@meta.data
df <- x1_meta %>%
  group_by(bin_x,xenium_clusters_p2) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_x) %>%
  mutate(prop = n / sum(n))

fn=paste0(output_PATH, target_tissue,"_",label,"_cell_cluster_dist_x_",id[i],".pdf")
p=ggplot(df,aes(x=bin_x,y=prop,fill=xenium_clusters_p2))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_brewer(palette="Set2")
pdf(fn)
print(p)
dev.off()
############

x1@meta.data$bin_y <- cut(
  x1@meta.data$y_coor,
  breaks = 100,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE
)

x1_meta=x1@meta.data
df <- x1_meta %>%
  group_by(bin_y,xenium_clusters_p2) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


#ggplot(df,aes(x=bin,y=prop,fill=xenium_clusters_p2))+geom_bar(stat = "identity") +
#  ylab("Proportion") +
#  theme_minimal()

fn=paste0(output_PATH, target_tissue,"_",label,"_cell_cluster_dist_y_",id[i],".pdf")
p=ggplot(df,aes(x=bin_y,y=prop,fill=xenium_clusters_p2))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_brewer(palette="Set2")
pdf(fn)
print(p)
dev.off()

#########


#for(i in 1:length(id)){

xenium_obj=xenium_data_list[[i]]

#meta=x1@meta.data #merged_obj@meta.data[merged_obj@meta.data$slide_id==id[i],]
xenium_obj@meta.data$xenium_clusters_p2="other"
xenium_obj@meta.data[x1_meta$barcode,]$xenium_clusters_p2=x1_meta$xenium_clusters_p2

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",id[i],".pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "xenium_clusters_p2", border.size = 0, cols = "polychrome", axes=TRUE)
#+ theme(
#  axis.line = element_line(),
#  axis.text = element_text(),
#  axis.title = element_text()
#) #,axis=TRUE)
print(p)
dev.off()
}




#}
#summary(lm(df[df$xenium_clusters_p2==0,]$prop~df[df$xenium_clusters_p2==0,]$bin))
