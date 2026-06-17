

library(Seurat)
library(dplyr)
library(polychrome)
####outdir <- "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/harmony"
args <- commandArgs(trailingOnly = TRUE)

#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/neuron_xenium.combined_annot.rds")
#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/nonneuron_xenium.combined_annot.rds")
#parts <- "nonneuron"


od="20250620_PP_comb" #args[1]
target_tissue="ONONH" #args[2]
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="3" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
#target_tissue="ONONH"
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0)

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", ".rds"))


merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_1"] <- "PP_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_2"] <- "PP_2"


merged_obj@meta.data$slide="snRNA"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_1"]="xenium_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_2"]="xenium_2"


table(merged_obj@meta.data[,c("majorclass","harmony_clusters_p3")])

tb=table(merged_obj@meta.data[,c("harmony_clusters_p3","majorclass")])
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_heatmap.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(p3_clu=rownames(tb),ct=max_colnames)

merged_obj@meta.data$harmony_anno="Unk"

for( i in anno$p3_clu){
merged_obj@meta.data$harmony_anno[merged_obj@meta.data$harmony_clusters_p3==i]=anno[anno$p3_clu==i,]$ct
}


xenium_data_list=readRDS(paste0(output_PATH, target_tissue, "_xenium_data_list_cutoff20.rds"))

for(i in 1:length(xenium_data_list)){

xenium_obj=xenium_data_list[[i]]
#xenium_obj@meta.data$harmony_anno[xenium_obj@meta.data$slide_id==paste0("PP_",i)]
#xenium_obj@meta.data$harmony_anno=merged_obj@meta.data$harmony_anno[(colnames(xenium_obj))&&(merged_obj@meta.data$slide_id==paste0("PP_",i))]

meta=merged_obj@meta.data[merged_obj@meta.data$sampleid==paste0("PP_",i),]
#xenium_obj@meta.data$harmony_anno=meta$harmony_anno[(colnames(xenium_obj))]
xenium_obj@meta.data$harmony_anno=meta[xenium_obj$barcode,"harmony_anno"]

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"
 
pdf(paste0(output_PATH, target_tissue, "_xenium_spatial_RCTD_annot_subclass_PP_",i,".pdf"), width = 50, height = 35)
ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0, cols = "polychrome",axis=TRUE)
print(p)
#ImageDimPlot(object = xenium_obj, group.by = "deconv_annot", border.size = 0, cols = "polychrome")
#ImageDimPlot(object = xenium_obj, group.by = "spot_class", border.size = 0, cols = "polychrome")
dev.off()
}
