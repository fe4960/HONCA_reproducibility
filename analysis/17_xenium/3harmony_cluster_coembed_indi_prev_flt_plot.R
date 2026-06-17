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


#library(harmony)


label="astro"
cut=20
sam="BCM_22_0698_ONH_RNA"
sc="subclass2"
xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,"_",sam,"_",sc,".rds"))


xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))


id=unique(xenium.combined@meta.data$slide_id)

library(ggplot2)
#############for(i in 2:length(id)){
############n=i-1
##########xenium_obj=xenium_data_list[[n]]
#########slide_id=unique(xenium_data_list[[n]]@meta.data$slide_id)
######merged_obj1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]
#######xenium_obj@meta.data$harmony_anno="other"
#######xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$harmony_anno=merged_obj1@meta.data$harmony_anno
########DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"
########pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_",id[i], "_",label, "_rpca_flt",cut,"_",sam,"_",sc,".pdf"), width = 50, height = 35)
########p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0,   cols = "Set2", axes=TRUE, dark.background=TRUE )
########print(p)
##########dev.off()
##########}

#####


x1_meta_all=NULL

for(i in 2:length(id)){


#x1@meta.data$x_coor=x1@images$fov$centroids$x

n=i-1
xenium_obj=xenium_data_list[[n]]
slide_id=unique(xenium_data_list[[n]]@meta.data$slide_id)

#x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]
x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

xenium_obj=xenium_obj[,x1@meta.data$barcode]
xenium_obj@meta.data$harmony_anno=x1@meta.data$harmony_anno

#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==4,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==5,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==6,]$xenium_clusters_p2=0

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

#pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",id[i],"_only_set3.pdf"), width = 50, height = 35)
pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_",id[i], "_",label, "_rpca_flt",cut,"_",sam,"_",sc,".pdf"), width = 50, height = 35)
#p=ImageDimPlot(object = xenium_obj, group.by = "xenium_clusters_p2", border.size = 0, cols = "polychrome", axes=TRUE)
p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0, cols = "Set3", axes=TRUE, dark.background=FALSE)
print(p)
dev.off()



x1@meta.data$y_coor=x1@images$fov$centroids$y


if((n==2) || (n==5)){
x1@meta.data$y_coor=-x1@images$fov$centroids$y
}

#library(ggplot2)
############



x1@meta.data$bin_y <- cut(
  x1@meta.data$y_coor,
  breaks = 50,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE

)

x1_meta=x1@meta.data

x1_meta_all=rbind(x1_meta_all,x1_meta)

}


#x1_meta_all[x1_meta_all$xenium_clusters_p2==4,]$xenium_clusters_p2=3
#x1_meta_all[x1_meta_all$xenium_clusters_p2==5,]$xenium_clusters_p2=4
#x1_meta_all[x1_meta_all$xenium_clusters_p2==6,]$xenium_clusters_p2=0

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

fn=paste0(output_PATH, target_tissue,"_",label,cut,"_",sam,"_",sc,"_cell_cluster_dist_y_comb.pdf")
p=ggplot(df,aes(x=bin_y,y=prop,fill=harmony_anno))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_brewer(palette="Set3")
pdf(fn)
print(p)
dev.off()

  

#========

