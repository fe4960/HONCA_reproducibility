library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)
#library(arrow)
library(RColorBrewer)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

cut=as.numeric(args[2])

########
#####process xenium data

label=args[1]

fc=as.numeric(args[3])
clu=args[4]
sc=args[5]
idx=args[6]
label1=args[7]
xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_6_xenium_data_list_comb_cutoff20.rds"))


xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))


id=unique(xenium.combined@meta.data$slide_id)
id=na.omit(id)
x1_meta_all=NULL

for(i in 1:length(id)){

xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]



if(idx == "n"){
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id, "_RCTD_results_nonneuron.rds"))
RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc,"_",label1, "_RCTD_results_nonneuron.rds"))
}else{
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron",idx,".rds"))
RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc,"_",label1, "_RCTD_results_nonneuron", idx, ".rds"))
}
x1@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
x1@meta.data[,"spot_class"]=RCTD@results$results_df[x1@meta.data$barcode,]$spot_class
x1=x1[,x1@meta.data[,"rctd1"]!="NA"]
x1=x1[,x1@meta.data$spot_class=="singlet"]

xenium_obj=xenium_obj[,x1@meta.data$barcode]
xenium_obj@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
xenium_obj@meta.data[,"spot_class"]=RCTD@results$results_df[x1@meta.data$barcode,]$spot_class

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

##########pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_",sc,"_",idx,"rctd1_only.pdf"), width = 50, height = 35)

###########p=ImageDimPlot(object = xenium_obj, group.by = "rctd1", border.size = 1,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "polychrome")
#p=ImageDimPlot(object = xenium_obj, group.by = "rctd1", border.size = 1,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "Set1") #oligo

##########print(p)
#########dev.off()



pdf(paste0(output_PATH, target_tissue, "_",label,"_cellsize_rpca_flt",cut,"_fc_",fc,"_",sc,"_slide_",slide_id,"_",idx,"_", clu,"_", label1,".pdf"), width = 15, height = 5)
#FeaturePlot(x1, reduction = "umap", features="cell_size",split.by="slide")
print(VlnPlot(x1, features="cell_size", group.by = "rctd1", pt.size =0))
print(VlnPlot(x1, features="nCount_Xenium", group.by = "rctd1", pt.size =0))
print(VlnPlot(x1, features="nFeature_Xenium", group.by = "rctd1", pt.size =0))

dev.off()


x1@meta.data$y_coor=x1@images$fov$centroids$x


########if((i==2) || (i==5)){
#########x1@meta.data$y_coor=-x1@images$fov$centroids$x
#}

#library(ggplot2)
############



x1@meta.data$bin_y <- cut(
  x1@meta.data$y_coor,
  breaks = 70,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE
)

x1_meta=x1@meta.data

x1_meta_all=rbind(x1_meta_all,x1_meta)

}



library(rlang)

library(dplyr)



df <- x1_meta_all %>%
  group_by(bin_y, rctd1   ) %>%

#  group_by(bin_y, !!sym(clu)  ) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


fn=paste0(output_PATH, target_tissue,"_",label,"_cutoff",cut,"_fc_",fc,"_",clu,"_cell_cluster_dist_y_comb_",sc,"_",label1,"_",idx,"_rmSlide2_5.pdf")
paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))

p=ggplot(df,aes(x=bin_y, y=prop, fill=rctd1  ))+geom_bar(stat = "identity") +
ylab("Proportion") + theme_minimal(base_size = 25) + scale_fill_manual( values= paired_more) # c("purple1", "lightblue",  "yellow",  "pink" ))
#scale_fill_brewer(palette="Set3")
pdf(fn, width=9)
print(p)
dev.off()

#########


