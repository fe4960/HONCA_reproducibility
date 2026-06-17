library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)
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

cut=args[2]

########
#####process xenium data

label=args[1]

fc=args[3]
clu=args[4]
sam=args[5]
sc=args[6]
################
#od="20250620_PP_comb" #args[1]
#target_tissue="ONONH_comb" #args[2]
############
#path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
#wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
#setwd(wd)
#output_PATH <- paste0(wd, "/",target_tissue, "/")
#cut=50


xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))


xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))

#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

for(i in 1:length(id)){


xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)


RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id, "_RCTD_results_nonneuron.rds"))

x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

xenium_obj=xenium_obj[,x1@meta.data$barcode]

xenium_obj@meta.data[,clu]=x1@meta.data[,clu]
xenium_obj@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
xenium_obj@meta.data[,"rctd2"]=RCTD@results$results_df[x1@meta.data$barcode,]$second_type

#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==4,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==5,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==6,]$xenium_clusters_p2=0

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",slide_id,".rds"))
xenium_obj@meta.data[,"harmony_clusters_p6"]=merged_obj@meta.data[paste0(slide_id,"_",x1@meta.data$barcode),]$harmony_clusters_p6
xenium_obj@meta.data[,"harmony_anno"]=merged_obj@meta.data[paste0(slide_id,"_",x1@meta.data$barcode),]$harmony_anno


tb=table(xenium_obj@meta.data[,c("rctd1",clu)])
max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(clu=rownames(tb),ct=max_colnames)


tb=table(xenium_obj@meta.data[,c("harmony_clusters_p6",clu)])
max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(clu=rownames(tb),ct=max_colnames)

tb=table(xenium_obj@meta.data[,c("harmony_anno",clu)])

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"





pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu,"_only.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = clu, border.size = NA,  axes=TRUE, dark.background=TRUE, cols = "polychrome", border.color=NA)
print(p)
dev.off()


pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "rctd1", border.size = NA,  axes=TRUE, dark.background=TRUE,  border.color=NA)
print(p)
dev.off()

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_harmony_anno_only.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = NA,  axes=TRUE, dark.background=TRUE,  border.color=NA)
print(p)
dev.off()


}

