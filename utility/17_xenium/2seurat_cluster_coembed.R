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
###saveRDS(xenium.combined, paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

#ref.combined=readRDS("/dfs3b/ruic20_lab/junw42/human_ret_anc/data/snRNA/ONONH_snRNA_subclass_5k.rds")
#meta=read.table("/dfs3b/ruic20_lab/junw42/human_ret_anc/data/snRNA/ONONH_snRNA_subclass_5k.obs.gz",header=T,sep="\t",row.names=1)

#/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz

ref.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)


#df=table(meta[meta$age_year>=60,]$sampleid)
#names(df[df>=3000])
#sample=c('A230245_Atlanta_ONH_RNA','BCM_22_0389_ONH_RNA','BCM_23_0334_ONH_RNA','GSM7553434','GSM7553443','GSM7553444','MMD_19_D008_ONH_RNA','MMD_D009_13_ONH_RNA')
#sample=c('A230245_Atlanta_ONH_RNA','BCM_22_0200_ON_RNA_s2','BCM_22_0389_ONH_RNA','BCM_22_0784_ON_RNA','BCM_23_0329_ONH_RNA','BCM_23_0334_ONH_RNA','MMD_19_D008_ONH_RNA','MMD_D009_13_ONH_RNA')
ref.combined <- AddMetaData(ref.combined, metadata = meta)
#ref.combined<- subset(x = ref.combined, subset = sampleid %in% sample)

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

ref.anchors <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, reduction = "rpca")
ref.combined <- IntegrateData(anchorset = ref.anchors)

DefaultAssay(ref.combined) <- "integrated"
ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30)
ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:30)
ref.combined <- FindClusters(ref.combined, resolution = 0.5)

saveRDS(ref.combined, paste0(output_PATH, target_tissue,"_ref.combined_5k_plus_panel.rds"))

pdf(paste0(output_PATH, target_tissue,"_ref.combined_5k_plus_panel.pdf"), width = 10, height = 5)
DimPlot(ref.combined)
DimPlot(ref.combined, group.by = "majorclass", label = T)
DimPlot(ref.combined, group.by = "sampleid", label = T)
dev.off()




###########

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(ref.combined)
DefaultAssay(xenium.combined)

obj_list <- list("ref" = ref.combined, "xenium" = xenium.combined)

features <- rownames(xenium.combined)
features <- intersect(features, rownames(ref.combined))

coembed.anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca")
coembed.combined <- IntegrateData(anchorset = coembed.anchors, new.assay.name = "coembed")
DefaultAssay(coembed.combined) <- "coembed"

coembed.combined <- ScaleData(coembed.combined, verbose = FALSE)
coembed.combined <- RunPCA(coembed.combined, npcs = 30, verbose = FALSE)
coembed.combined <- RunUMAP(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindNeighbors(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindClusters(coembed.combined, resolution = 0.5, cluster.name = "p5_clusters")


saveRDS(coembed.combined, paste0(output_PATH, target_tissue, "_coembed.combined_inte_ref_xenium.rds"))
coembed.combined@meta.data$stat=Idents(coembed.combined)
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_ref_xenium.pdf"), width = 6, height = 5)
DimPlot(coembed.combined, label = T, group.by = "p5_clusters") + ggtitle("Xenium Clusters p5")
DimPlot(coembed.combined, label = T, group.by = "ident") + ggtitle("ref or query")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium.pdf"), width = 12, height = 5)
DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
dev.off()

