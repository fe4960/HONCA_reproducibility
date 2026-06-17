library(future)
plan("multicore")  # or "multisession" if you're on Windows
options(future.globals.maxSize = 80 * 1024^3)  # 50 GB

library(Seurat)
library(ggplot2)
library(dplyr)
library(FNN)
library(harmony)
source("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/from_Tingting/ispatial/iSpatial.R")


args <- commandArgs(trailingOnly = TRUE)

# read xenium data
#slide_name <- "Round2_Slide04"
slide_name = "ONH_1" #<- args[1]
#print(slide_name)
outdir <- paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/imputation/",slide_name, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]

wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")


integrated_merFISH=readRDS(paste0(outdir, "/", slide_name, "_integrated_merFISH_metaxy.rds"))

#library(Seurat)

integrated_merFISH@images <- list()
coords <- integrated_merFISH@meta.data[, c("x","y")]
coords <- coords[colnames(integrated_merFISH), ]

centroids <- CreateCentroids(coords = coords)
fov <- CreateFOV(centroids)

# 4. attach to object
integrated_merFISH@images <- list(fov = fov)

#fov <- CreateFOV(coords = as.matrix(coords), type = "centroids")
#integrated_merFISH[["fov"]] <- fov


get_gene_xenium <- function(xenium.combined, gene, type){
  expr_mat <- xenium.combined@assays$enhanced$data[rownames(xenium.combined@assays$enhanced$data) %in% gene, xenium.combined@meta.data$harmony_anno == type]
  return(mean(expr_mat))
}

marker_list=c("GAD1","SLC4A11","GRM6","ARR3","PTPRB","BICC1","ONECUT1","CD69","F13A1", "MLANA","RGR","CD74","NOTCH3", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ")

#marker_list <- c("Atf3", "Calca", "Adra2a", "Bmpr1b",
#                 "Oprk1", "Smr2", "Sstr2", "Mrgpra3", "Mrgprb4", "Trpv1", "Mrgprd",
#                 "Ntrk3", "S100a16", "Ntrk2", "Sst","Th",
#                 "Trpm8")
#types <- sort(unique(integrated_merFISH@meta.data$harmony_anno))
types=c("AC","Astrocyte","BC","Cone","Endothelial_cell","Fibroblast","HC","Immune_cell","Macrophage","Melanocyte","MG","Microglia","Mural_cell","Oligodendrocyte","Oligodendrocyte_precursor_cell","RGC","RPE","Rod", "Schwann_cell")
mean_list <- list()
for (i in c(1:length(marker_list))){
  gene <- marker_list[i]
  mean_list[[i]] <- sapply(types, function(x){get_gene_xenium(integrated_merFISH, gene, x)})
}
hm_df <- as.data.frame(do.call(cbind, mean_list))
colnames(hm_df) <- marker_list
rownames(hm_df) <- types
saveRDS(hm_df, paste0(outdir, "/", slide_name, "_hm_df.rds"))
library(viridis)
library(pheatmap)
pdf(paste0(outdir, "/", slide_name, "_imputed_marker_heatmap.pdf"), width = 8, height = 6)
p <- pheatmap(hm_df, scale = "column", cluster_cols = F, cluster_rows = F, color = viridis(n = 10, alpha = 1,
                                                                                           begin = 0, end = 1, option = "viridis"), border_color = "black")
print(p)
dev.off()

#############



integrated_merFISH <- FindVariableFeatures(integrated_merFISH, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated_merFISH)
integrated_merFISH <- ScaleData(integrated_merFISH, features = all.genes)
integrated_merFISH <- RunPCA(integrated_merFISH, npcs = 30, features = VariableFeatures(integrated_merFISH))
integrated_merFISH <- RunUMAP(integrated_merFISH, reduction = "pca", dims = 1:30)
integrated_merFISH <- FindNeighbors(integrated_merFISH, dims = 1:30, reduction = "pca")
integrated_merFISH <- FindClusters(integrated_merFISH, resolution = 1, cluster.name = "imputed_w1_by_slide")

#saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH_clustering.rds"))
saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH_clustering_metaxy.rds"))

pdf(paste0(outdir, "/", slide_name, "_umap.pdf"), width = 8, height = 5)
DimPlot(integrated_merFISH)
DimPlot(integrated_merFISH, group.by = "harmony_anno", label = T, repel = T)
dev.off()
#pdf(paste0(outdir, "/", slide_name, "_neuron_umap_tglabel.pdf"), width = 10, height = 5)
#DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T, split.by = "tglabel")
#dev.off()
#pdf(paste0(outdir, "/", slide_name, "_neuron_umap_branchlabel.pdf"), width = 15, height = 5)
#DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T, split.by = "branchlabel")
#dev.off()



