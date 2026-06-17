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

label=args[1]

cut=args[2]
fc=args[3]
sc=args[4]
idx=args[5]
#####process xenium data


#clu="xenium_clusters_p3" #args[4]
sam="BCM_22_0047_ON_ONH_RNA_ds" #args[5]


###########################xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))
#########################xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_8_xenium_data_list_comb_cutoff20.rds"))


xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_8_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))


#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

#ii=c(6,8)
#for(i in 1:length(id)){
for(i in 6:7){



xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)



x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

xenium_obj=xenium_obj[,x1@meta.data$barcode]

fn=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron",idx,".rds")
if(idx == "n"){
fn=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds")
}
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))

RCTD=readRDS(fn)

xenium_obj@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
xenium_obj@meta.data[,"rctd2"]=RCTD@results$results_df[x1@meta.data$barcode,]$second_type
xenium_obj@meta.data[,"spot_class"]=RCTD@results$results_df[x1@meta.data$barcode,]$spot_class

xenium_obj=xenium_obj[,xenium_obj@meta.data[,"rctd1"]!="NA"]
xenium_obj=xenium_obj[,xenium_obj@meta.data$spot_class=="singlet"]

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"







pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_",sc,"_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only1_dark_color.pdf"), width = 22, height = 12)
#######p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = c("purple1",  "yellow" )) # oligo #, crop=TRUE, overlap=TRUE,  alpha=1) # oligo

######p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = c("purple1", "lightblue",  "yellow",  "pink", "green" )) #, crop=TRUE, overlap=TRUE,  alpha=1) #astro

#########p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "Paired") #, crop=TRUE, overlap=TRUE,  alpha=1)


p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = c("purple1", "lightblue",  "yellow",  "pink", "green" )) #, crop=TRUE, overlap=TRUE,  alpha=1)

#####p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = c("purple1", "lightblue",  "yellow",  "pink" )) #, crop=TRUE, overlap=TRUE,  alpha=1)
print(p)
dev.off()


#######p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "polychrome") #, crop=TRUE, overlap=TRUE,  alpha=1)



}

