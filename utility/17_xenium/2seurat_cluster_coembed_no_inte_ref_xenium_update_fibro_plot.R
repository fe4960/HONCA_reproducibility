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


library(harmony)




merged_obj=readRDS(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony.rds"))

merged_obj <- FindClusters(merged_obj, resolution = 1.2, cluster.name = "w12_clusters")

merged_obj <- FindClusters(merged_obj, resolution = 1.5, cluster.name = "w15_clusters")

merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "w2_clusters")

merged_obj1=subset(merged_obj, subset = subclass!="unknown") 

#merged_obj1@meta.data$subclass1=paste(strsplit(merged_obj1@meta.data$subclass, "_")[[1]][1:2], collapse="_")

get_first_two <- function(s) {
  words <- unlist(strsplit(s, "_"))
  paste(words[1:2], collapse = "_")
}

merged_obj1@meta.data$subclass1=sapply( merged_obj1@meta.data$subclass , get_first_two)

merged_obj@meta.data$subclass1=sapply( merged_obj@meta.data$subclass , get_first_two)


tb=table(merged_obj1@meta.data[,c("p8_clusters","subclass")])
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_heatmap.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

tb=table(merged_obj1@meta.data[,c("p8_clusters","subclass1")])
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_heatmap_subclass1.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

tb=table(merged_obj1@meta.data[,c("w1_clusters","subclass1")])
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_heatmap_subclass1_w1.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()


tb=table(merged_obj1@meta.data[,c("w12_clusters","subclass1")])
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_heatmap_subclass1_1.2.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

tb=table(merged_obj1@meta.data[,c("w15_clusters","subclass1")])
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_heatmap_subclass1_1.5.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

tb=table(merged_obj1@meta.data[,c("w2_clusters","subclass1")])
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_heatmap_subclass1_2.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()



max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(p1_2_clu=rownames(tb),ct=max_colnames)
#anno[anno$p1_2_clu==11,"ct"]="Fibro_perivascular"
merged_obj@meta.data$harmony_anno="Unk"

for( i in anno$p1_2_clu){
	
#merged_obj@meta.data$harmony_anno[merged_obj@meta.data$harmony_clusters_1_2==i]=anno[anno$p1_2_clu==i,]$ct
#merged_obj@meta.data$harmony_anno[merged_obj@meta.data$p8_clusters==i]=anno[anno$p1_2_clu==i,]$ct

merged_obj@meta.data$harmony_anno[merged_obj@meta.data$w2_clusters==i]=anno[anno$p1_2_clu==i,]$ct

}


merged_obj <- merged_obj %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

#pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_p8_umap_slide.pdf"), width = 20, height = 5)
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_2_umap_slide.pdf"), width = 20, height = 5)

DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno", split.by="slide_id", label=TRUE)+ NoLegend()
dev.off()


pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_2_umap_ident.pdf"), width = 10, height = 5)

#pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_p8_umap_ident.pdf"), width = 10, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno", split.by="ident", label =TRUE) + NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_2_umap_ident_subclass1.pdf"), width = 10, height = 5)

#pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_p8_umap_ident_subclass1.pdf"), width = 10, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="subclass1", split.by="ident", label =TRUE) + NoLegend()
DimPlot(merged_obj, reduction = "umap", group.by="subclass1", split.by="ident", label =FALSE) + NoLegend()

dev.off()




merged_obj2=merged_obj[,merged_obj@meta.data$nFeature_Xenium>=5]

#pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_p8_umap_anno_nFea5_p8.pdf"), width = 6, height = 5)
pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_2_umap_anno_nFea5_2.pdf"), width = 6, height = 5)


DimPlot(merged_obj2, reduction = "umap", group.by="harmony_anno", split.by="ident", label =TRUE) + NoLegend()

dev.off()


#saveRDS(merged_obj, file=paste0(output_PATH, target_tissue, "fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_p8", ".rds"))
saveRDS(merged_obj, file=paste0(output_PATH, target_tissue, "fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_anno_2", ".rds"))


cut=20


library(ggplot2)

#####slide=c("PP_1", "PP_2", "ONH_1", "ONH_2", "ONH_3")
######v=c("", ".2", ".3", ".4", ".5")

########for(i in 1:5){
######xenium_obj=subset(merged_obj, subset = slide_id == slide[i])
######DefaultBoundary(xenium_obj[[paste0("fov",v[i])]]) <-  "segmentation"
######pdf(paste0(output_PATH, target_tissue, "_fibro_annot_subclass1_", slide[i], "_p8.pdf"), width = 50, height = 35)
#######p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0,   cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )
########print(p)
#######p=ImageDimPlot(object = xenium_obj, group.by = "p8_clusters", border.size = 0,   cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )
#######print(p)
########dev.off()
#######3}

#####

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))

cut=20
xenium.combined=merged_obj   #readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1", ".rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

#id=unique(xenium.combined@meta.data$slide_id)

library(ggplot2)
for(i in 1:length(xenium_data_list)){

xenium_obj=xenium_data_list[[i]]

slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

merged_obj1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]
xenium_obj=xenium_obj[,merged_obj1@meta.data$barcode]


xenium_obj@meta.data$p8_clusters="other"

xenium_obj@meta.data$harmony_anno="other"
xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$harmony_anno=merged_obj1@meta.data$harmony_anno
xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$p8_clusters=merged_obj1@meta.data$p8_clusters

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

#pdf(paste0(output_PATH, target_tissue, "_fibro_annot_subclass1_", slide_id, "_p8.pdf"), width = 50, height = 35)
pdf(paste0(output_PATH, target_tissue, "_fibro_annot_subclass1_", slide_id, "_2.pdf"), width = 50, height = 35)

p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0,   cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )
print(p)

p=ImageDimPlot(object = xenium_obj, group.by = "p8_clusters", border.size = 0,   cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )
print(p)
dev.off()
}



#========

