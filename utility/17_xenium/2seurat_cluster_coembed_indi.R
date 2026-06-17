library(Seurat)
library(ggplot2)
library(dplyr)
#library(arrow)
args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
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

samplelist=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_sampleid_sum_top20_flt4",header=T,sep="\t")





# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0) 

snRNA_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/" 
ref.combined=readRDS(paste0(snRNA_dir,"Oligodendrocyte_subclass_seurat_simple.rds"))
meta=read.table(paste0(snRNA_dir,"Oligodendrocyte_subclass_seurat_simple.obs.gz"),header=T,sep=",",row.names=1)


ref.combined <- AddMetaData(ref.combined, metadata = meta)

ref.combined = subset(ref.combined, subset = sampleid %in% samplelist$sampleid)


xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))


features <- rownames(xenium.combined)
#########
Idents(ref.combined) <- "sampleid"
ref.list <- SplitObject(ref.combined, split.by = "sampleid")

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- intersect(features, rownames(ref.combined))
features <- c(features, SelectIntegrationFeatures(object.list = ref.list))
features <- unique(features)

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#ref.anchors <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, reduction = "cca")

ref.anchors <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, reduction = "rpca")
ref.combined <- IntegrateData(anchorset = ref.anchors)

DefaultAssay(ref.combined) <- "integrated"
ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30)
ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:30)
ref.combined <- FindClusters(ref.combined, resolution = 0.5)

#saveRDS(ref.combined, paste0(output_PATH, target_tissue,"_OLIGO_ref.combined_487_plus_7panel.rds"))
saveRDS(ref.combined, paste0(output_PATH, target_tissue,"_OLIGO_ref.combined_487_plus_4panel_rpca.rds"))
pdf(paste0(output_PATH, target_tissue,"_OLIGO_ref.combined_487_plus_4panel_rpca.pdf"), width = 10, height = 5)

#pdf(paste0(output_PATH, target_tissue,"_OLIGO_ref.combined_487_plus_7panel.pdf"), width = 10, height = 5)
DimPlot(ref.combined)
DimPlot(ref.combined, group.by = "subclass", label = T)
DimPlot(ref.combined, group.by = "sampleid", label = T)
dev.off()

#sl=c("BCM_22_0784_ON_RNA","BCM_22_0849_ON_RNA","MMD_19_D008_ON_RNA","MMD_22_16608_ON_RNA")
#ref.combined1 = subset(ref.combined, subset = !(sampleid %in% sl))

#DimPlot(ref.combined1, group.by = "sampleid", label = T)
