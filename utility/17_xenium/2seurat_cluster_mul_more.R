library(Seurat)
library(ggplot2)
library(dplyr)
#install.packages("arrow")
library(arrow)
args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="6" #args[6]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0) 
matching_dirs0 <- list.dirs(path, full.names = TRUE, recursive = FALSE)
xenium_data_list <- list()
matching_dirs=matching_dirs0[grepl(kw,matching_dirs0)]

###########


for (i in c(1:length(matching_dirs))){
  dir <- paste0(kw, i)
  xenium_path_full <- matching_dirs[i]
  print(xenium_path_full)
  # adjust to new version of xenium ranger, which dropping transcripts.csv.gz in the output folder
  #transcripts <- read_parquet(file.path(xenium_path_full, "transcripts.parquet"))
  #write.csv(transcripts, gzfile(file.path(xenium_path_full, "transcripts.csv.gz")), row.names = FALSE)
#####  xenium.obj <- LoadXenium(xenium_path_full, fov = "fov")
  ##########
  
  data <- ReadXenium(xenium_path_full, outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
names(data)
## continue the regular LoadXenium
  segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
	xmo <- CreateSeuratObject(
	  counts = data$matrix[["Gene Expression"]], 
	  assay = "Xenium")
	xmo[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
	xmo[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
	xmo[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
	xmo[["fov"]] <- coords

	print(xmo)
   xenium.obj = xmo

  ##########
  print(dim(xenium.obj))
  xenium.obj@meta.data$slide_id <- dir #unlist(strsplit(dir, split = "_"))[1]

  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
 md=basename(normalizePath(matching_dirs[i]))
    cell_id_PATH <- paste0(cb_dir,"/",md)
    cell_list_list <- list()
    j <- 1
      cell_list_file <- paste0(cell_id_PATH, "/",  "Selection_",id,"_cells_stats.csv")
      cell_list <- read.table(cell_list_file, header = T,comment.char = "#",sep=",")
      cell_list_list[[j]] <- cell_list #[,1]
    cell_list <- do.call(rbind, cell_list_list)

  xenium.obj@meta.data$barcode <- rownames(xenium.obj@meta.data)
  xenium.obj <- subset(xenium.obj, subset = barcode %in% cell_list[,1])
  print(dim(xenium.obj))
  
  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
  
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 20)
  print(dim(xenium.obj))
   print("start") 
    xenium.obj@meta.data$cell_size <- sapply(rownames(xenium.obj@meta.data), FUN = function(x){xenium.obj@images$fov$segmentation[x]@polygons[[1]]@area})
   print("end")
#  xenium.obj@meta.data$cell_size <- sapply(rownames(xenium.obj@meta.data), FUN = function(x){xenium.obj@images$fov$segmentation[x]@polygons[[1]]@area})
  
  xenium_data_list[[i]] <- xenium.obj
  names(xenium_data_list)[i] <- dir #unlist(strsplit(dir, split = "_"))[1]
  
}

############################
path0="20250711__194041__HumAd_ONH-PeriRetRPECho_20250711" #args[3]
kw="ONH_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="6" #args[6]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0) 
matching_dirs0_1 <- list.dirs(path, full.names = TRUE, recursive = FALSE)
###to sort directory name alphabetically:
matching_dirs0_1=sort(matching_dirs0_1)
matching_dirs_1=matching_dirs0_1[grepl(kw,matching_dirs0_1)]
xenium_data_list_1 <- list()


for (i in c(1:length(matching_dirs_1))){

#for (i in c(1:length(matching_dirs_1))){
  dir <- paste0(kw, i)
  xenium_path_full <- matching_dirs_1[i]
  print(xenium_path_full)
  # adjust to new version of xenium ranger, which dropping transcripts.csv.gz in the output folder
  transcripts <- read_parquet(file.path(xenium_path_full, "transcripts.parquet"))
  write.csv(transcripts, gzfile(file.path(xenium_path_full, "transcripts.csv.gz")), row.names = FALSE)
  ######xenium.obj <- LoadXenium(xenium_path_full, fov = "fov")

  data <- ReadXenium(xenium_path_full, outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
names(data)
## continue the regular LoadXenium
  segmentations.data <- list(
  centroids = CreateCentroids(data$centroids),
  segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(
  coords = segmentations.data, 
  type = c("segmentation", "centroids"), 
  molecules = data$microns, 
  assay = "Xenium")
	xmo <- CreateSeuratObject(
	  counts = data$matrix[["Gene Expression"]], 
	  assay = "Xenium")
	xmo[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
	xmo[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
	xmo[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
	xmo[["fov"]] <- coords

	print(xmo)
   xenium.obj = xmo

 
#########

  print(dim(xenium.obj))
  xenium.obj@meta.data$slide_id <- dir #unlist(strsplit(dir, split = "_"))[1]

  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
 md=basename(normalizePath(matching_dirs_1[i]))
    cell_id_PATH <- paste0(cb_dir,"/",md)
    cell_list_list <- list()
    j <- 1
      cell_list_file <- paste0(cell_id_PATH, "/",  "Selection_",id,"_cells_stats.csv")
      cell_list <- read.table(cell_list_file, header = T,comment.char = "#",sep=",")
      cell_list_list[[j]] <- cell_list #[,1]
    cell_list <- do.call(rbind, cell_list_list)

  xenium.obj@meta.data$barcode <- rownames(xenium.obj@meta.data)
  xenium.obj <- subset(xenium.obj, subset = barcode %in% cell_list[,1])
  print(dim(xenium.obj))
  
  print(quantile(xenium.obj$nCount_Xenium, seq(0, 1, 0.1)))
  print(quantile(xenium.obj$nFeature_Xenium, seq(0, 1, 0.1)))
  
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 20)
  print(dim(xenium.obj))
  
  xenium.obj@meta.data$cell_size <- sapply(rownames(xenium.obj@meta.data), FUN = function(x){xenium.obj@images$fov$segmentation[x]@polygons[[1]]@area})
  
  xenium_data_list_1[[i]] <- xenium.obj
  names(xenium_data_list_1)[i] <- dir #unlist(strsplit(dir, split = "_"))[1]
  
}

xenium_data_list_comb=c(xenium_data_list, xenium_data_list_1)

saveRDS(xenium_data_list_comb, paste0(output_PATH, target_tissue, "_6_xenium_data_list_comb_cutoff20.rds"))


xenium_data_list_comb <- lapply(X = xenium_data_list_comb, FUN = function(x) {
  x <- NormalizeData(x)
})

features <- rownames(xenium.obj)
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
xenium.combined <- FindClusters(xenium.combined, resolution = 0.8, cluster.name = "xenium_clusters_p8")
xenium.combined <- FindClusters(xenium.combined, resolution = 1, cluster.name = "xenium_clusters_1")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.2, cluster.name = "xenium_clusters_p2")

dim(xenium.combined) # >20: TM 818; CB 36134; surface: 28505 ONONH: 478 49290
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_6_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_6_xenium_umap_p2_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T) + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3") + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, group.by = "slide_id") + ggtitle("Sections")
FeaturePlot(xenium.combined, label = T, features = "nCount_Xenium", max.cutoff = "q90")
FeaturePlot(xenium.combined, label = T, features = "nFeature_Xenium")
FeaturePlot(xenium.combined, label = T, features = "cell_size")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_6_xenium_umap_1_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_1") + ggtitle("Xenium Clusters 1")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p5") + ggtitle("Xenium Clusters p5")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p8") + ggtitle("Xenium Clusters p8")

dev.off()

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_6_markers_p2.csv"))

table(xenium.combined@meta.data$xenium_clusters_p3)

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_6_markers_p3.csv"))


markers <- FindAllMarkers(xenium.combined, group.by = "xenium_clusters_p5", only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_6_markers_p5.csv"))


#marker_genes <- c("MLANA", "SCN7A", "PECAM1", "CCL21", "PTPRC", "ACTA2", "DES", "MME", "POU6F2", "KRT5", "PAX6", "RGS5", "PDGFRB")
marker_genes=c("GFAP","AC092957.1","MECOM","VWF","BICC1","GABRG3","HKDC1","F13A1","PAX3","KCNQ3","ADAM28","MYH11","CARMN","RNF220","CTNNA3","NXPH1","CSMD1","BCL11B","PYHIN1","CD247","KIT","SKAP1","BANK1")
pdf(paste0(output_PATH, target_tissue, "_6_xenium_vln_dot_p5_cutoff20.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = marker_genes, pt.size = 0)
DotPlot(xenium.combined, features = marker_genes)
dev.off()
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_6_xenium.combined_cutoff20.rds"))

opc=c("CSMD1","TNR","OLIG1","OLIG2","OPCML","SNTG1")

peri=c("RGS5","ABCC9","TRPC4")

oligo=c("MOG","MOBP","ST18")

mk=c(opc,peri, oligo)

xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_6_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_6_xenium_vln_dot_p2_cutoff20.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = mk, pt.size = 0)
DotPlot(xenium.combined, features = mk)
dev.off()

