library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
#library(arrow)
args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)

ref.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.rds")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple_sort.obs.gz",header=T,sep=",",row.names=1)



ref.combined <- AddMetaData(ref.combined, metadata = meta)

sub="t"

sample=c("BCM_22_0784_ON_RNA","MMD_23_17738_ON_RNA","BCM_23_0491_ONH_RNA","MMD_19_D008_ONH_RNA","MMD_23_22486_ONH_RNA","BCM_22_0890_ON_RNA")

if(sub == "t"){
	ref.combined=subset(ref.combined, subset = sampleid %in% sample)
}

xenium.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_fibro_xenium.combined_cutoff20_fc_5.rds")
features <- rownames(xenium.combined)
#########
Idents(ref.combined) <- "sampleid"

features <- intersect(features, rownames(ref.combined))

ref.combined=ref.combined[features,]
ref.combined=NormalizeData(ref.combined)

xenium.combined@meta.data$subclass1="unknown"
xenium.combined@meta.data$subclass="unknown"


ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"

DefaultAssay(ref.combined)
DefaultAssay(xenium.combined)


obj_list <- list("ref" = ref.combined, "xenium" = xenium.combined)
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


coembed.combined <- FindClusters(coembed.combined, resolution = 0.6, cluster.name = "p6_clusters")
coembed.combined <- FindClusters(coembed.combined, resolution = 0.7, cluster.name = "p7_clusters")
coembed.combined <- FindClusters(coembed.combined, resolution = 0.8, cluster.name = "p8_clusters")

coembed.combined <- RunUMAP(coembed.combined, reduction = "harmony", dims = 1:30)

saveRDS(coembed.combined, paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_subset.rds"))

pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_cluster_subset.pdf"), width = 6, height = 5)


DimPlot(coembed.combined, label = T, group.by = "p8_clusters") + ggtitle("Xenium Clusters p8")

DimPlot(coembed.combined, label = T, group.by = "p7_clusters") + ggtitle("Xenium Clusters p7")

DimPlot(coembed.combined, label = T, group.by = "p6_clusters") + ggtitle("Xenium Clusters p6")

DimPlot(coembed.combined, label = T, group.by = "p5_clusters") + ggtitle("Xenium Clusters p5")
DimPlot(coembed.combined, label = T, group.by = "p3_clusters") + ggtitle("Xenium Clusters p3")
DimPlot(coembed.combined, label = T, group.by = "p2_clusters") + ggtitle("Xenium Clusters p2")
DimPlot(coembed.combined, label = T, group.by = "w1_clusters") + ggtitle("Xenium Clusters 1")

DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
dev.off()


pdf(paste0(output_PATH, target_tissue, "_fibro_ref_xenium_no_inte_cutoff20_fc_5_harmony_cluster_facet_subset.pdf"), width = 12, height = 5)
DimPlot(coembed.combined, group.by = "subclass", reduction = "umap",  label=T, split.by="ident")
dev.off()


