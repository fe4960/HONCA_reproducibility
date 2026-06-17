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

ref.combined=readRDS(args[1]) #readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table(args[2], header=T,sep=",",row.names=1) #read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)


#df=table(meta[meta$age_year>=60,]$sampleid)
#names(df[df>=3000])
#sample=c('A230245_Atlanta_ONH_RNA','BCM_22_0389_ONH_RNA','BCM_23_0334_ONH_RNA','GSM7553434','GSM7553443','GSM7553444','MMD_19_D008_ONH_RNA','MMD_D009_13_ONH_RNA')
#sample=c('A230245_Atlanta_ONH_RNA','BCM_22_0200_ON_RNA_s2','BCM_22_0389_ONH_RNA','BCM_22_0784_ON_RNA','BCM_23_0329_ONH_RNA','BCM_23_0334_ONH_RNA','MMD_19_D008_ONH_RNA','MMD_D009_13_ONH_RNA')
ref.combined <- AddMetaData(ref.combined, metadata = meta)

ref.combined = subset(ref.combined, subset = sampleid == args[5])
#ref.combined<- subset(x = ref.combined, subset = sampleid %in% sample)
#ref.combined <- subset(x = ref.combined, downsample = as.numeric(args[5]))

#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_cutoff20.rds"))
label=args[4]

for(i in length(xenium_data_list)){

bc=read.table(file=paste0(args[3],"_",i),header=T,row.names=1)
xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)
}




features <- rownames(xenium_data_list[[1]])
features <- intersect(features, rownames(ref.combined))

xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")
xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)

DefaultAssay(xenium.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#xenium.combined <- ScaleData(xenium.combined, verbose = FALSE)
#xenium.combined <- RunPCA(xenium.combined, npcs = 30, verbose = FALSE)
#xenium.combined <- RunUMAP(xenium.combined, reduction = "pca", dims = 1:30)

#features <- rownames(xenium.combined)
#########
# Idents(ref.combined) <- "sampleid"
# split the dataset into a list of two seurat objects (stim and CTRL)
######ref.list <- SplitObject(ref.combined, split.by = "sampleid")

# normalize and identify variable features for each dataset independently
######ref.list <- lapply(X = ref.list, FUN = function(x) {
#######  x <- NormalizeData(x)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
######})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
#features <- rownames(xenium.combined)
#######features <- c(features, SelectIntegrationFeatures(object.list = ref.list))
#######features <- unique(features)

########ref.list <- lapply(X = ref.list, FUN = function(x) {
#######  x <- ScaleData(x, features = features, verbose = FALSE)
#######  x <- RunPCA(x, features = features, verbose = FALSE)
#########})

#######3ref.anchors <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, reduction = "rpca")
########ref.combined <- IntegrateData(anchorset = ref.anchors)

#######3DefaultAssay(ref.combined) <- "integrated"
######3ref.combined <- ScaleData(ref.combined, verbose = FALSE)
######ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
#######ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30)
#######ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:30)
########ref.combined <- FindClusters(ref.combined, resolution = 0.5)

#######saveRDS(ref.combined, paste0(output_PATH, target_tissue,"_ref.combined_5k_plus_panel.rds"))

#######pdf(paste0(output_PATH, target_tissue,"_ref.combined_5k_plus_panel.pdf"), width = 10, height = 5)
#########DimPlot(ref.combined)
#######DimPlot(ref.combined, group.by = "majorclass", label = T)
#######3DimPlot(ref.combined, group.by = "sampleid", label = T)
########dev.off()

ref.combined=ref.combined[features,]
ref.combined=NormalizeData(ref.combined)


###########
# co-embedding
#xenium.combined <- readRDS("DRG_Xenium/xenium.combined.rds")
#drg_ref.combined <- readRDS("DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(ref.combined)
DefaultAssay(xenium.combined)

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
coembed.combined <- ScaleData(coembed.combined, verbose = FALSE)
coembed.combined <- RunPCA(coembed.combined, npcs = 30, verbose = FALSE)
coembed.combined <- RunUMAP(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindNeighbors(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindClusters(coembed.combined, resolution = 0.5, cluster.name = "p5_clusters")
#coembed.combined <- FindClusters(coembed.combined, resolution = 0.8, cluster.name = "Rp8_clusters")


saveRDS(coembed.combined, paste0(output_PATH, target_tissue, "_coembed.combined_inte_ref_xenium_no_inte_",label,".rds"))
coembed.combined@meta.data$stat=Idents(coembed.combined)
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_ref_xenium_no_inte_",label,".pdf"), width = 6, height = 5)
DimPlot(coembed.combined, label = T, group.by = "p5_clusters") + ggtitle("Xenium Clusters p5")
DimPlot(coembed.combined, label = T, group.by = "ident") + ggtitle("ref or query")
#DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_",label,".pdf"), width = 10, height = 5)
DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
dev.off()

#library(SeuratDisk)
#filename=paste0(output_PATH, target_tissue, "_coembed.combined_inte.h5Seurat")
#SaveH5Seurat(coembed.combined, filename = filename)
#Convert(filename, dest = "h5ad")

coembed.combined=readRDS(paste0(output_PATH, target_tissue, "_coembed.combined_inte_ref_xenium_no_inte_",label,".rds"))
#pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_ct.pdf"), width = 12, height = 5)
#DimPlot(coembed.combined, reduction = "umap", split.by = "ident", group.by="majorclass",label=T)
#dev.off()

#coembed.combined@meta.data$[coembed.combined@meta.data$ident=="Xenium",]$sampleid="Unknown"

#pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_ct.pdf"), width = 12, height = 5)
#DimPlot(coembed.combined, reduction = "umap", split.by = "ident", group.by="majorclass",label=T)+ NoLegend()
#dev.off()

#coembed.combined@meta.data$[coembed.combined@meta.data$ident=="Xenium",]$sampleid="Unknown"


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_ct_",label,".pdf"), width = 5, height = 5)
DimPlot(coembed.combined, reduction = "umap", group.by="subclass",label=T)+ NoLegend()
dev.off()


#coembed.combined@meta.data[coembed.combined@meta.data$slide_id=="PP_1",]$sampleid="PP_1"
#coembed.combined@meta.data[coembed.combined@meta.data$slide_id=="PP_2",]$sampleid="PP_2"

coembed.combined@meta.data$sampleid[coembed.combined@meta.data$slide_id == "PP_1"] <- "PP_1"
coembed.combined@meta.data$sampleid[coembed.combined@meta.data$slide_id == "PP_2"] <- "PP_2"

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_sampleid_",label,".pdf"), width = 5, height = 5)
DimPlot(coembed.combined, reduction = "umap", group.by="sampleid")+ NoLegend()
dev.off()

coembed.combined@meta.data$slide="snRNA"
coembed.combined@meta.data$slide[coembed.combined@meta.data$slide_id == "PP_1"]="xenium_1"
coembed.combined@meta.data$slide[coembed.combined@meta.data$slide_id == "PP_2"]="xenium_2"

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_slide_",label,".pdf"), width = 15, height = 5)
DimPlot(coembed.combined, reduction = "umap", split.by="slide")+ NoLegend()
dev.off()

