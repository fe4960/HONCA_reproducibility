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

#cut=50

########
#####process xenium data

#label="oligo"

################
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
#cut=50


xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))


obs=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/sn_st_scvi_harmony_obs.txt.gz", sep="\t", header=T, row.names=1)

score=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/sn_st_scvi_harmony_score.txt.gz", sep="\t", header=T, row.names=1)
row_max <- apply(score, 1, max)

obs=obs[row_max >=0.7,]

x1_meta_all=NULL

for(i in 1:5){


xenium_obj=xenium_data_list[[i]]


slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

obs1=obs[obs$sampleid==slide_id,]

xenium_obj=xenium_obj[,obs1$barcode]

xenium_obj@meta.data$lr_celltype="other"
xenium_obj@meta.data[obs1$barcode,]$lr_celltype=obs1$lr_celltype
xenium_obj=xenium_obj[,!(xenium_obj@meta.data$lr_celltype %in% c("Pigmented_cell","Adipocyte"))]

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"
x1=xenium_obj

x1@meta.data$y_coor=x1@images$fov$centroids$y

if((i==2) || (i==5)){
x1@meta.data$y_coor=-x1@images$fov$centroids$y
}

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
  group_by(bin_y, lr_celltype) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


library(RColorBrewer)

paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))

fn=paste0(output_PATH, target_tissue,"_cell_cluster_dist_y_comb_scvi.pdf")
p=ggplot(df,aes(x=bin_y,y=prop,fill=lr_celltype))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_manual(values=paired_more)
pdf(fn)
print(p)
dev.off()


