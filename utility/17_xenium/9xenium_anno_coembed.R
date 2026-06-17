library(Seurat)
library(ggplot2)
library(dplyr)
#library(arrow)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" 
target_tissue="ONONH_comb" 
############
path0="20250620__203138__Human_70y_ASRetina_20250620" 
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/"
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

cut=50
label="oligo"


xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,".rds"))


merged_obj@meta.data$xenium_cluster_p2="NA"

merged_obj@meta.data[colnames(xenium.combined),]$xenium_cluster_p2=xenium.combined@meta.data$xenium_clusters_p2

table(merged_obj@meta.data[,c("subclass","xenium_cluster_p2")])

tb=table(merged_obj@meta.data[,c("xenium_cluster_p2","subclass")])
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_heatmap_xenium_",label,"_rpca_flt",cut,".pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(xenium_cluster_p2=rownames(tb),ct=max_colnames)

merged_obj@meta.data$harmony_anno_xenium="Unk"

for( i in anno$xenium_cluster_p2){
merged_obj@meta.data$harmony_anno_xenium[merged_obj@meta.data$xenium_cluster_p2==i]=anno[anno$xenium_cluster_p2==i,]$ct
}


merged_obj@meta.data$slide="snRNA"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_1"]="xenium_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_2"]="xenium_2"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_1"] <- "xenium_3"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_2"] <- "xenium_4"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_3"] <- "xenium_5"

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_anno_rpca_flt",cut,"_xenium.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno",label=F,split.by="slide") #+ NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_subclass_rpca_flt",cut,"_xenium.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="subclass",label=F,split.by="slide")
DimPlot(merged_obj, reduction = "umap", group.by="xenium_cluster_p2",label=F,split.by="slide") 
#+ NoLegend()
dev.off()

