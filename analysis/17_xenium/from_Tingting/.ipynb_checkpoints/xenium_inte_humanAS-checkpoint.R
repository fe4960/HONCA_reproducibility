# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(arrow)

getwd()
setwd("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/human_eye_xenium_batch250225/")

target_tissue <- "Surface"
output_PATH <- paste0("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/human_eye_xenium_batch250225/basic_analysis/", target_tissue, "/")

# read xenium data
path <- "/dfs3b/ruic20_lab/rawdata/Xenium/20250224_Human_Anterior_Segment_panel_G33B7Z/20250220__211510__HuAdultAnteriorSegment_20250220/"
matching_dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)

xenium_data_list <- list()
for (i in c(1:length(matching_dirs))){
  dir <- paste0("sec0", i)
  xenium_path_full <- matching_dirs[i]
  print(xenium_path_full)
  
  # adjust to new version of xenium ranger, which dropping transcripts.csv.gz in the output folder
  #transcripts <- read_parquet(file.path(xenium_path_full, "transcripts.parquet"))
  #write.csv(transcripts, gzfile(file.path(xenium_path_full, "transcripts.csv.gz")), row.names = FALSE)
  
  xenium.obj <- LoadXenium(xenium_path_full, fov = "fov")
  print(dim(xenium.obj))
  xenium.obj@meta.data$slide_id <- dir #unlist(strsplit(dir, split = "_"))[1]
  
  #pdf(paste0("basic_analysis/", dir, "_basic.pdf"))
  #p <- ImageFeaturePlot(xenium.obj, fov = "fov", features = "nCount_Xenium", max.cutoff = "q90") +
  #  theme(
  #    panel.grid = element_blank(),
  #    panel.border = element_blank(),
  #    axis.ticks = element_blank(),
  #    axis.text = element_blank()
  #  )
  #print(p)
  #p <- ImageFeaturePlot(xenium.obj, fov = "fov", features = "nFeature_Xenium", max.cutoff = "q90") +
  #  theme(
  #    panel.grid = element_blank(),
  #    panel.border = element_blank(),
  #    axis.ticks = element_blank(),
  #    axis.text = element_blank()
  #  )
  #print(p)
  #dev.off()
  
  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
  
  if (target_tissue == "Surface"){
    cell_id_PATH <- paste0("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/human_eye_xenium_batch250225/basic_analysis/")
    cell_list_list <- list()
    j <- 1
    for (tissue in c("Cornea", "Limbus", "Sclera")){
      cell_list_file <- paste0(cell_id_PATH, tissue, "/",  "sec0", i, "_", tissue, "_cells_id.txt")
      cell_list <- read.csv(cell_list_file, header = FALSE)
      cell_list_list[[j]] <- cell_list
      j <- j + 1
    }
    cell_list <- do.call(rbind, cell_list_list)

  } else {
    cell_id_PATH <- paste0("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/human_eye_xenium_batch250225/basic_analysis/", target_tissue, "/")
    cell_list_file <- paste0(cell_id_PATH, "sec0", i, "_", target_tissue, "_cells_id.txt")
    cell_list <- read.csv(cell_list_file, header = FALSE)
  }
  

  xenium.obj@meta.data$barcode <- rownames(xenium.obj@meta.data)
  xenium.obj <- subset(xenium.obj, subset = barcode %in% cell_list$V1)
  print(dim(xenium.obj))
  
  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
  
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 20)
  print(dim(xenium.obj))
  
  xenium.obj@meta.data$cell_size <- sapply(rownames(xenium.obj@meta.data), FUN = function(x){xenium.obj@images$fov$segmentation[x]@polygons[[1]]@area})
  
  xenium_data_list[[i]] <- xenium.obj
  names(xenium_data_list)[i] <- unlist(strsplit(dir, split = "_"))[1]
  
}

saveRDS(xenium_data_list, paste0(output_PATH, target_tissue, "_xenium_data_list_cutoff20.rds"))


xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- NormalizeData(x)
})

features <- rownames(xenium.obj)
xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")
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

dim(xenium.combined) # >20: TM 818; CB 36134; surface: 28505
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p1_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T) + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3") + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, group.by = "slide_id") + ggtitle("Sections")
FeaturePlot(xenium.combined, label = T, features = "nCount_Xenium", max.cutoff = "q90")
FeaturePlot(xenium.combined, label = T, features = "nFeature_Xenium")
FeaturePlot(xenium.combined, label = T, features = "cell_size")
dev.off()

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_markers_p3.csv"))

table(xenium.combined@meta.data$xenium_clusters_p3)

marker_genes <- c("MLANA", "SCN7A", "PECAM1", "CCL21", "PTPRC", "ACTA2", "DES", "MME", "POU6F2", "KRT5", "PAX6", "RGS5", "PDGFRB")
VlnPlot(xenium.combined, features = marker_genes, pt.size = 0)
DotPlot(xenium.combined, features = marker_genes)

xenium.combined@meta.data$xenium_clusters_p3_celltype <- recode(xenium.combined@meta.data$xenium_clusters_p3, 
                                                             "0"= "Fibroblasts", "1"= "Epithelium", "2"= "Epithelium", "3"= "Endothelium/Pericytes", "4"= "Epithelium",
                                                             "5"= "Macrophages", "6"= "Corneal Stromal Fibroblasts", "7"= "TM/CB Fibroblasts", "8"= "Melanocytes",
                                                             "9"= "Corneal Endothelium", "10"= "NK/T Cells", "11" = "Schwann Cells", "12" = "Lymphatic Endothelium",
                                                             "13" = "Mast Cells", "14" = "Smooth Muscle Cells")
xenium.combined@meta.data$xenium_clusters_p3_majorclass <- recode(xenium.combined@meta.data$xenium_clusters_p3, 
                                                                "0"= "Fibroblasts", "1"= "Epithelium", "2"= "Epithelium", "3"= "Endothelium/Pericytes", "4"= "Epithelium",
                                                                "5"= "Immune Cells", "6"= "Fibroblasts", "7"= "Fibroblasts", "8"= "Melanocytes",
                                                                "9"= "Corneal Endothelium", "10"= "Immune Cells", "11" = "Schwann Cells", "12" = "Lymphatic Endothelium",
                                                                "13" = "Immune Cells", "14" = "Smooth Muscle Cells")
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

xenium.combined <- readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p1_cutoff20_annot.pdf"), width = 10, height = 5)
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3_celltype", repel = T) + ggtitle("Xenium Annotation Using Markers")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3_majorclass", repel = T) + ggtitle("Xenium Annotation Using Markers")
dev.off()

marker_genes <- c("MLANA", "SCN7A", "PECAM1", "CCL21", "PTPRC", "C1QA", "CD96", "CPA3", "ACTA2", "DES", "LUM", "MME", "PI16", "POU6F2", "KRT5", "PAX6", "RGS5")
VlnPlot(xenium.combined, features = marker_genes, pt.size = 0, group.by = "xenium_clusters_p3_celltype", stack = T)
DotPlot(xenium.combined, features = marker_genes, group.by = "xenium_clusters_p3_celltype", scale = T, scale.by = "radius")
VlnPlot(xenium.combined, features = "SLC6A6", pt.size = 0, group.by = "xenium_clusters_p3_celltype")

pdf(paste0(output_PATH, target_tissue, "_xenium_cutoff20_annot_vln.pdf"), width = 15, height = 5)
VlnPlot(xenium.combined, features = marker_genes, pt.size = 0, group.by = "xenium_clusters_p3_celltype", stack = T)
dev.off()

write.csv(xenium.combined@meta.data, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20_annot_obs.csv"))
saveRDS(xenium.combined@meta.data, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20_annot_obs.rds"))

table(xenium.combined$xenium_clusters_p3_celltype)/nrow(xenium.combined@meta.data)*100
table(xenium.combined$xenium_clusters_p3_majorclass)/nrow(xenium.combined@meta.data)*100

epi_xenium <- subset(xenium.combined, subset = xenium_clusters_p3_celltype == "Epithelium")
VlnPlot(epi_xenium, features = c("KRT5", "KRT3", "KRT15", "KRT4"), pt.size = 0, group.by = "xenium_clusters_1", stack = T)
VlnPlot(epi_xenium, features = c("KRT3", "KRT15"), pt.size = 0, group.by = "xenium_clusters_p3", stack = T)
VlnPlot(epi_xenium, features = c("KRT3", "KRT15"), pt.size = 0, group.by = "xenium_clusters_p5", stack = T)
VlnPlot(epi_xenium, features = c("KRT3", "KRT15"), pt.size = 0, group.by = "xenium_clusters_1", stack = T)


#================================================================================================
# reference
raw_ref <- readRDS("basic_analysis/Surface/Surface_BCM_22_0500.rds")
table(raw_ref@meta.data$majorclass)

# subset to panel genes
mtx <- raw_ref@assays$RNA$counts[features,]
dim(mtx)
ref_panel <- CreateSeuratObject(mtx)
ref_panel@meta.data <- raw_ref@meta.data
ref_panel <- NormalizeData(ref_panel)

Idents(ref_panel) <- "majorclass"
markers <- FindAllMarkers(ref_panel, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_ref_markers.csv"))
top_10_vector <- markers %>%
  group_by(cluster) %>%  # Group by columnC
  #arrange(desc(columnB)) %>%  # Arrange in descending order of columnB
  slice_head(n = 10) %>%  # Select the top 10 rows for each group
  pull(gene)

top_10_vector <- unique(top_10_vector)
pdf(paste0(output_PATH, target_tissue, "_DEG_dotplot_ref_panel_majorclass.pdf"), width = 30, height = 8)
DotPlot(ref_panel, features = top_10_vector) + RotatedAxis() + xlab("Top10 DEGs") + ylab("")
dev.off()

Idents(ref_panel) <- "celltype"
markers <- FindAllMarkers(ref_panel, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_ref_markers_celltype.csv"))
top_10_vector <- markers %>%
  group_by(cluster) %>%  # Group by columnC
  #arrange(desc(columnB)) %>%  # Arrange in descending order of columnB
  slice_head(n = 5) %>%  # Select the top 10 rows for each group
  pull(gene)

top_10_vector <- unique(top_10_vector)
pdf(paste0(output_PATH, target_tissue, "_DEG_dotplot_ref_panel_celltype.pdf"), width = 30, height = 8)
DotPlot(ref_panel, features = top_10_vector) + RotatedAxis() + xlab("Top5 DEGs") + ylab("")
dev.off()

ref_panel <- NormalizeData(ref_panel)
#ref_panel <- FindVariableFeatures(ref_panel)
ref_panel <- ScaleData(ref_panel)
ref_panel <- RunPCA(ref_panel, assay = "RNA", npcs = 30, features = features)
ref_panel <- FindNeighbors(ref_panel, assay = "RNA", reduction = "pca", dims = 1:30)
ref_panel <- FindClusters(ref_panel, cluster.name = "seurat_cluster", resolution = 0.3)
ref_panel <- RunUMAP(ref_panel, reduction = "pca", reduction.name = "umap", return.model = T, dims = 1:30)

pdf(paste0(output_PATH, target_tissue, "_ref_panel_umaps_uninte.pdf"))
DimPlot(ref_panel, reduction = "umap", label = F) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref_panel, reduction = "umap", label = T) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref_panel, reduction = "umap", label = F, group.by = "majorclass") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref_panel, reduction = "umap", label = T, group.by = "majorclass") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref_panel, reduction = "umap", label = T, group.by = "celltype", repel = T) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref_panel, reduction = "umap", label = F, group.by = "sampleid") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
#DimPlot(ref_panel, reduction = "umap", label = F, group.by = "Sample_info") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
dev.off()

Idents(ref_panel) <- "sampleid"
ref_list <- SplitObject(ref_panel, split.by = "sampleid")

# normalize and identify variable features for each dataset independently
ref_list <- lapply(X = ref_list, FUN = function(x) {
  x <- NormalizeData(x)
  #x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

ref_list <- lapply(X = ref_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

ref.anchors <- FindIntegrationAnchors(object.list = ref_list, anchor.features = features, reduction = "rpca")
ref.combined <- IntegrateData(anchorset = ref.anchors)

DefaultAssay(ref.combined) <- "integrated"
ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30)
ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:30)
ref.combined <- FindClusters(ref.combined, resolution = 0.5)

saveRDS(ref.combined, paste0(output_PATH, target_tissue, "_ref.combined_inte_panel_genes.rds"))

pdf(paste0(output_PATH, target_tissue, "_ref_panel_umaps.pdf"), width = 8, height = 5)
DimPlot(ref.combined, reduction = "umap", label = F) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref.combined, reduction = "umap", label = T) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref.combined, reduction = "umap", label = F, group.by = "majorclass") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref.combined, reduction = "umap", label = T, group.by = "majorclass") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref.combined, reduction = "umap", label = T, group.by = "celltype", repel = T) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref.combined, reduction = "umap", label = F, group.by = "sampleid") + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_ref_panel_umaps_celltype.pdf"), width = 15, height = 8)
DimPlot(ref.combined, reduction = "umap", label = F, group.by = "celltype", repel = T) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
DimPlot(ref.combined, reduction = "umap", label = T, group.by = "celltype", repel = T) + ggtitle(paste0(target_tissue, " - Xenium Panel Set")) + labs(subtitle = "BCM_22_0500")
dev.off()

table(ref.combined$majorclass)/nrow(ref.combined@meta.data)*100
table(ref.combined$celltype)/nrow(ref.combined@meta.data)*100





