

library(Seurat)
library(dplyr)
library(ggplot2)
#library(polychrome)
####outdir <- "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/harmony"
args <- commandArgs(trailingOnly = TRUE)

#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/neuron_xenium.combined_annot.rds")
#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/nonneuron_xenium.combined_annot.rds")
#parts <- "nonneuron"



options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="5" #args[6]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

cut=50

label="oligo"

merged_obj=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))


id=unique(merged_obj@meta.data$slide_id)

#xenium_data_list=readRDS(paste0(output_PATH, target_tissue, "_xenium_data_list_cutoff20.rds"))

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))

for(i in 1:length(id)){

xenium_obj=xenium_data_list[[i]]

meta=merged_obj@meta.data[merged_obj@meta.data$slide_id==id[i],]
xenium_obj@meta.data$xenium_clusters_p2="other"
xenium_obj@meta.data[meta$barcode,]$xenium_clusters_p2=meta$xenium_clusters_p2

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"
 
pdf(paste0(output_PATH, target_tissue, "_xenium_spatial_",label,"_annot_subclass_",id[i],".pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "xenium_clusters_p2", border.size = 0, cols = "polychrome")+ theme(
  axis.line = element_line(),
  axis.text = element_text(),
  axis.title = element_text()
) #,axis=TRUE)
print(p)
dev.off()
}
