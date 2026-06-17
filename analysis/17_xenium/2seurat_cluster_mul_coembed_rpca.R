library(Seurat)
library(ggplot2)
library(dplyr)
library(arrow)

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
#dir.create(output_PATH, recursive = TRUE)
# read xenium data


#xenium_data_list_comb=readRDS(paste0(output_PATH, target_tissue, "_8_xenium_data_list_comb_cutoff20.rds"))

xenium_data_list_comb=readRDS(paste0(output_PATH, target_tissue, "_8_xenium_data_list_comb_cutoff20.rds"))


xenium_data_list_comb <- lapply(X = xenium_data_list_comb, FUN = function(x) {
#  x <- NormalizeData(x)
colnames(x)=paste0(x$slide_id,"_",x$barcode)
x <- NormalizeData(x)
#return(x)
})


#xenium_data_list_comb <- lapply(X = xenium_data_list_comb, FUN = function(x) {
#  x <- NormalizeData(x)
#})

features <- rownames(xenium_data_list_comb[[1]])
xenium_data_list_comb <- lapply(X = xenium_data_list_comb, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list_comb, anchor.features = features, reduction = "rpca")
xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)

DefaultAssay(xenium.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
xenium.combined <- ScaleData(xenium.combined, verbose = FALSE)
xenium.combined <- RunPCA(xenium.combined, npcs = 30, verbose = FALSE)
xenium.combined <- RunUMAP(xenium.combined, reduction = "pca", dims = 1:30)
xenium.combined <- FindNeighbors(xenium.combined, reduction = "pca", dims = 1:30)
xenium.combined <- FindClusters(xenium.combined, resolution = 0.3, cluster.name = "xenium_clusters_p3")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.5, cluster.name = "xenium_clusters_p5")
xenium.combined <- FindClusters(xenium.combined, resolution = 1, cluster.name = "xenium_clusters_1")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.2, cluster.name = "xenium_clusters_p2")

dim(xenium.combined) # >20: TM 818; CB 36134; surface: 28505 ONONH: 478 49290
#saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_8_xenium.combined_cutoff20.rds"))


saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_8_xenium.combined_cutoff20.rds"))


pdf(paste0(output_PATH, target_tissue, "_8_xenium_umap_p2_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T) + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3") + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, group.by = "slide_id") + ggtitle("Sections")
FeaturePlot(xenium.combined, label = T, features = "nCount_Xenium", max.cutoff = "q90")
FeaturePlot(xenium.combined, label = T, features = "nFeature_Xenium")
FeaturePlot(xenium.combined, label = T, features = "cell_size")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_8_xenium_umap_1_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_1") + ggtitle("Xenium Clusters 1")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p5") + ggtitle("Xenium Clusters p5")

dev.off()

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_8_markers_p2.csv"))

table(xenium.combined@meta.data$xenium_clusters_p3)

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_8_markers_p3.csv"))


markers <- FindAllMarkers(xenium.combined, group.by = "xenium_clusters_p5", only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_8_markers_p5.csv"))


#marker_genes <- c("MLANA", "SCN7A", "PECAM1", "CCL21", "PTPRC", "ACTA2", "DES", "MME", "POU6F2", "KRT5", "PAX6", "RGS5", "PDGFRB")
marker_genes=c("GFAP","AC092957.1","MECOM","VWF","BICC1","GABRG3","HKDC1","F13A1","PAX3","KCNQ3","ADAM28","MYH11","CARMN","RNF220","CTNNA3","NXPH1","CSMD1","BCL11B","PYHIN1","CD247","KIT","SKAP1","BANK1")
pdf(paste0(output_PATH, target_tissue, "_8_xenium_vln_dot_p5_cutoff20.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = marker_genes, pt.size = 0)
DotPlot(xenium.combined, features = marker_genes)
dev.off()
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_8_xenium.combined_cutoff20.rds"))

opc=c("CSMD1","TNR","OLIG1","OLIG2","OPCML","SNTG1")

peri=c("RGS5","ABCC9","TRPC4")

oligo=c("MOG","MOBP","ST18")

mk=c(opc,peri, oligo)

####xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_8_xenium_vln_dot_p2_cutoff20.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = mk, pt.size = 0)
DotPlot(xenium.combined, features = mk)
dev.off()

