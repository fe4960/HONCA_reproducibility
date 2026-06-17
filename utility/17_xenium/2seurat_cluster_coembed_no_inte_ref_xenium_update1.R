library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

library(future)
#plan(sequential)
options(future.globals.maxSize = 50 * 1024^3)


#library(arrow)
args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
#path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0) 

ref.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)


ref.combined <- AddMetaData(ref.combined, metadata = meta)

#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))
xenium.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_8_xenium.combined_cutoff20.rds")
#xenium_data_list=readRDS(paste0(output_PATH, target_tissue, "_xenium_data_list_cutoff20.rds"))
xenium_data_list=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_8_xenium_data_list_comb_cutoff20.rds")
features <- rownames(xenium.combined)
#########
Idents(ref.combined) <- "sampleid"

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
#features <- rownames(xenium.combined)
features <- intersect(features, rownames(ref.combined))

ref.combined=ref.combined[features,]
ref.combined=NormalizeData(ref.combined)


#xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
#  x=x[features,]			   
#  x <- NormalizeData(x)
#  x@meta.data$sampleid=x@meta.data$slide_id
#  x@meta.data$majorclass="unknown"

#})



#Idents(xenium_data_list) <- "sampleid"
###########
# co-embedding
#xenium.combined <- readRDS("DRG_Xenium/xenium.combined.rds")
#drg_ref.combined <- readRDS("DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
#xenium_data_list@meta.data$ident <- "Xenium"
#Idents(xenium_data_list) <- "sampleid"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
#Idents(xenium_data_list) <- "ident"

DefaultAssay(ref.combined)
DefaultAssay(xenium.combined)
#DefaultAssay(xenium_data_list)


#obj_list <- list("ref" = ref.combined, "xenium" = xenium_data_list)

obj_list <- list("ref" = ref.combined, "xenium" = xenium.combined)
#obj_list <- lapply(X = obj_list, FUN = function(x) {
#  x <- NormalizeData(x)
#})
#features <- rownames(xenium.combined)
#features <- intersect(features, rownames(ref.combined))
obj_list <- lapply(X = obj_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

coembed.anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca")
coembed.combined <- IntegrateData(anchorset = coembed.anchors, new.assay.name = "coembed")
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(coembed.combined) <- "coembed"

# Run the standard workflow for visualization and clustering
coembed.combined <- ScaleData(coembed.combined, verbose = FALSE, vars.to.regress = "ident")
coembed.combined <- RunPCA(coembed.combined, npcs = 30, verbose = FALSE)


coembed.combined <- RunHarmony(coembed.combined, "ident")
coembed.combined <- FindNeighbors(coembed.combined, reduction = "harmony", dims = 1:30)
coembed.combined <- FindClusters(coembed.combined, resolution = 0.2, cluster.name = "p2_clusters")
coembed.combined <- FindClusters(coembed.combined, resolution = 0.3, cluster.name = "p3_clusters")

coembed.combined <- FindClusters(coembed.combined, resolution = 0.5, cluster.name = "p5_clusters")
coembed.combined <- FindClusters(coembed.combined, resolution = 1, cluster.name = "w1_clusters")

coembed.combined <- RunUMAP(coembed.combined, reduction = "harmony", dims = 1:30)
#coembed.combined <- FindClusters(coembed.combined, resolution = 0.8, cluster.name = "Rp8_clusters")


saveRDS(coembed.combined, paste0(output_PATH, target_tissue, "_8_coembed.combined_cutoff20_ref_xenium_no_inte_harmony.rds"))
#coembed.combined@meta.data$stat=Idents(coembed.combined)

pdf(paste0(output_PATH, target_tissue, "_8_xenium_umap_cutoff20_ref_xenium_no_inte.pdf"), width = 6, height = 5)
DimPlot(coembed.combined, label = T, group.by = "p5_clusters") + ggtitle("Xenium Clusters p5")
DimPlot(coembed.combined, label = T, group.by = "p3_clusters") + ggtitle("Xenium Clusters p3")
DimPlot(coembed.combined, label = T, group.by = "p2_clusters") + ggtitle("Xenium Clusters p2")
DimPlot(coembed.combined, label = T, group.by = "w1_clusters") + ggtitle("Xenium Clusters 1")

#DimPlot(coembed.combined, label = T, group.by = "ident") + ggtitle("ref or query")
DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
#DimPlot(coembed.combined, label = T, group.by = "majorclass") + ggtitle("ref or query")
#DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_8_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte.pdf"), width = 12, height = 5)
DimPlot(coembed.combined, group.by = "majorclass", reduction = "umap",  label=T)
dev.off()


