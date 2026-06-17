library(Seurat)
library(ggplot2)
library(dplyr)
#library(arrow)
args <- commandArgs(trailingOnly = TRUE)
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
matching_dirs0 <- list.dirs(path, full.names = TRUE, recursive = FALSE)
xenium_data_list <- list()
matching_dirs=matching_dirs0[grepl(kw,matching_dirs0)]
#md=matching_dirs

oligo1=read.table(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_1"),row.names=1,header=T)
oligo2=read.table(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_2"),row.names=1,header=T)

opc1=read.table(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc_1"),row.names=1,header=T)
opc2=read.table(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc_2"),row.names=1,header=T)

fibro1=read.table(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro_1"),row.names=1,header=T)
fibro2=read.table(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro_2"),row.names=1,header=T)


xenium_data_list=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_data_list_cutoff20.rds")



pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_1.pdf"), width=10, height=10)

ImageDimPlot(xenium_data_list[[1]], fov = "fov", molecules=c("CERCAM","TMEM144","ST18","CTNNA3"), axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",coord.fixed = FALSE,  nmols = 10000, cells=oligo1$x)
dev.off()
for (i in c(1:length(matching_dirs))){
  
  dir <- paste0(kw, i)
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
 md=basename(normalizePath(matching_dirs[i]))
#  if (target_tissue == "Surface"){
    cell_id_PATH <- paste0(cb_dir,"/",md)

 #    cell_id_PATH <- paste0(cb_dir,"/",matching_dirs[i])
    cell_list_list <- list()
    j <- 1
 #   for (tissue in c("Cornea", "Limbus", "Sclera")){
#      cell_list_file <- paste0(cell_id_PATH, tissue, "/",  "sec0", i, "_", tissue, "_cells_id.txt")
      cell_list_file <- paste0(cell_id_PATH, "/",  "Selection_",id,"_cells_stats.csv")

      cell_list <- read.table(cell_list_file, header = T,comment.char = "#",sep=",")
      cell_list_list[[j]] <- cell_list #[,1]
#      j <- j + 1
#    }
    cell_list <- do.call(rbind, cell_list_list)

#  } 
#else {
#    cell_id_PATH <- paste0("/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/human_eye_xenium_batch250225/basic_analysis/", target_tissue, "/")
#    cell_list_file <- paste0(cell_id_PATH, "sec0", i, "_", target_tissue, "_cells_id.txt")
#    cell_list <- read.csv(cell_list_file, header = FALSE)
#  }

xenium.obj@meta.data$barcode <- rownames(xenium.obj@meta.data)
  xenium.obj <- subset(xenium.obj, subset = barcode %in% cell_list[,1])
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

dim(xenium.combined) # >20: TM 818; CB 36134; surface: 28505 ONONH: 478 49290
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p2_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T) + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3") + ggtitle("Xenium Clusters")
DimPlot(xenium.combined, group.by = "slide_id") + ggtitle("Sections")
FeaturePlot(xenium.combined, label = T, features = "nCount_Xenium", max.cutoff = "q90")
FeaturePlot(xenium.combined, label = T, features = "nFeature_Xenium")
FeaturePlot(xenium.combined, label = T, features = "cell_size")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_1_cutoff20.pdf"), width = 6, height = 5)
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_1") + ggtitle("Xenium Clusters 1")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p5") + ggtitle("Xenium Clusters p5")

dev.off()

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_markers_p2.csv"))

table(xenium.combined@meta.data$xenium_clusters_p3)

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_markers_p3.csv"))


markers <- FindAllMarkers(xenium.combined, group.by = "xenium_clusters_p5", only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue, "_markers_p5.csv"))


#marker_genes <- c("MLANA", "SCN7A", "PECAM1", "CCL21", "PTPRC", "ACTA2", "DES", "MME", "POU6F2", "KRT5", "PAX6", "RGS5", "PDGFRB")
marker_genes=c("GFAP","AC092957.1","MECOM","VWF","BICC1","GABRG3","HKDC1","F13A1","PAX3","KCNQ3","ADAM28","MYH11","CARMN","RNF220","CTNNA3","NXPH1","CSMD1","BCL11B","PYHIN1","CD247","KIT","SKAP1","BANK1")
pdf(paste0(output_PATH, target_tissue, "_xenium_vln_dot_p5_cutoff20.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = marker_genes, pt.size = 0)
DotPlot(xenium.combined, features = marker_genes)
dev.off()
saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

opc=c("CSMD1","TNR","OLIG1","OLIG2","OPCML","SNTG1")

peri=c("RGS5","ABCC9","TRPC4")

oligo=c("MOG","MOBP","ST18")

mk=c(opc,peri, oligo)

xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_vln_dot_p2_cutoff20.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = mk, pt.size = 0)
DotPlot(xenium.combined, features = mk)
dev.off()

