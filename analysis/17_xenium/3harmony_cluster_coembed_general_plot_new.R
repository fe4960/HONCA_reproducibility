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




merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", "_fea.rds"))


#merged_obj <- FindClusters(merged_obj, resolution = 1.2, cluster.name = "harmony_clusters_1_2")
#merged_obj <- FindClusters(merged_obj, resolution = 1.5, cluster.name = "harmony_clusters_1_5")
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "harmony_clusters_2")

merged_obj <- merged_obj %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

table(merged_obj@meta.data[,c("majorclass","harmony_clusters_1_2")])
#table(merged_obj@meta.data[,c("majorclass","harmony_clusters_1_5")])
#table(merged_obj@meta.data[,c("majorclass","harmony_clusters_2")])

tb=table(merged_obj@meta.data[,c("harmony_clusters_1_2","majorclass")])
#tb=table(merged_obj@meta.data[,c("harmony_clusters_1_5","majorclass")])
#tb=table(merged_obj@meta.data[,c("harmony_clusters_2","majorclass")])

#pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_clu_2_ct_heatmap.pdf"), width = 7, height = 7)

#pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_clu_1_5_ct_heatmap.pdf"), width = 7, height = 7)
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_clu_1_2_ct_heatmap_fea.pdf"), width = 7, height = 7)

heatmap(tb)
dev.off()

max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(p1_2_clu=rownames(tb),ct=max_colnames)

merged_obj@meta.data$harmony_anno="Unk"

for( i in anno$p1_2_clu){
merged_obj@meta.data$harmony_anno[merged_obj@meta.data$harmony_clusters_1_2==i]=anno[anno$p1_2_clu==i,]$ct
}

merged_obj@meta.data[merged_obj@meta.data$harmony_anno=="T_cell",]$harmony_anno="Immune_cell"
merged_obj@meta.data[merged_obj@meta.data$harmony_anno=="Oligodendrocyte_precursor_cell",]$harmony_anno="OPC"


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_slide_ct_anno_fea.pdf"), width = 20, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno", split.by="slide", label=TRUE)+ NoLegend()
dev.off()



pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ident_ct_anno_fea.pdf"), width = 10, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno", split.by="ident", label =TRUE) + NoLegend()
dev.off()


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_clu_1_2_ct_anno_fea.pdf"), width = 10, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_clusters_1_2", split.by="ident", label =TRUE) + NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_majorclass_ct_anno_fea.pdf"), width = 10, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="majorclass", split.by="ident", label =TRUE) + NoLegend()
dev.off()



saveRDS(merged_obj, file=paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", "_fea.rds"))

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))

cut=20
xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", "_fea.rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)

library(ggplot2)
for(i in 1:(length(id)-1)){


xenium_obj=xenium_data_list[[i]]

slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

merged_obj1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]


#meta=x1@meta.data #merged_obj@meta.data[merged_obj@meta.data$slide_id==id[i],]

xenium_obj@meta.data$harmony_anno="other"
xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$harmony_anno=merged_obj1@meta.data$harmony_anno

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_",id[i],"_fea.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0,  alpha=0.8, cols = "polychrome",border.color=NA, axes=TRUE, dark.background=TRUE )

#+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(p)
dev.off()

}

#####

#========

