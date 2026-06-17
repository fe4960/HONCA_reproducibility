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

ref.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)


ref.combined <- AddMetaData(ref.combined, metadata = meta)


#####generate a integrated object for xenium data, the default assay of the integrated object should be "integrated" 
xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))


features <- rownames(xenium.combined)

#########only do Normalization for reference snRNA-seq, no integration across sampleid is needed. This approach may be more suitable for the reference from a small number of samples without too much batch effect.

#Idents(ref.combined) <- "sampleid"

features <- intersect(features, rownames(ref.combined))

ref.combined=ref.combined[features,]
ref.combined=NormalizeData(ref.combined)


###########
# co-embedding

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

DefaultAssay(coembed.combined) <- "coembed"

coembed.combined <- ScaleData(coembed.combined, verbose = FALSE)
coembed.combined <- RunPCA(coembed.combined, npcs = 30, verbose = FALSE)
coembed.combined <- RunUMAP(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindNeighbors(coembed.combined, reduction = "pca", dims = 1:30)
coembed.combined <- FindClusters(coembed.combined, resolution = 0.5, cluster.name = "p5_clusters")


saveRDS(coembed.combined, paste0(output_PATH, target_tissue, "_coembed.combined_inte_ref_xenium_no_inte.rds"))
coembed.combined@meta.data$stat=Idents(coembed.combined)
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_ref_xenium_no_inte.pdf"), width = 6, height = 5)
DimPlot(coembed.combined, label = T, group.by = "p5_clusters") + ggtitle("Xenium Clusters p5")
DimPlot(coembed.combined, label = T, group.by = "ident") + ggtitle("ref or query")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte.pdf"), width = 10, height = 5)
DimPlot(coembed.combined, reduction = "umap", split.by = "ident")
dev.off()

#####the following are for downstream analysis
coembed.combined=readRDS(paste0(output_PATH, target_tissue, "_coembed.combined_inte_ref_xenium_no_inte.rds"))


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_ct.pdf"), width = 5, height = 5)
DimPlot(coembed.combined, reduction = "umap", group.by="majorclass",label=T)+ NoLegend()
dev.off()



coembed.combined@meta.data$sampleid[coembed.combined@meta.data$slide_id == "PP_1"] <- "PP_1"
coembed.combined@meta.data$sampleid[coembed.combined@meta.data$slide_id == "PP_2"] <- "PP_2"

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_sampleid.pdf"), width = 5, height = 5)
DimPlot(coembed.combined, reduction = "umap", group.by="sampleid")+ NoLegend()
dev.off()

coembed.combined@meta.data$slide="snRNA"
coembed.combined@meta.data$slide[coembed.combined@meta.data$slide_id == "PP_1"]="xenium_1"
coembed.combined@meta.data$slide[coembed.combined@meta.data$slide_id == "PP_2"]="xenium_2"

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_slide.pdf"), width = 15, height = 5)
DimPlot(coembed.combined, reduction = "umap", split.by="slide")+ NoLegend()
dev.off()


oligo1=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_1")&(coembed.combined@meta.data$p5_clusters %in% c(1,12,15,18,23,28,29))]

oligo2=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_2")&(coembed.combined@meta.data$p5_clusters %in% c(1,12,15,18,23,28,29))]


opc1=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_1")&(coembed.combined@meta.data$p5_clusters==11)]
opc2=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_2")&(coembed.combined@meta.data$p5_clusters==11)]

fibro1=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_1")&(coembed.combined@meta.data$p5_clusters %in% c(4,7,2,17))]
fibro2=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_2")&(coembed.combined@meta.data$p5_clusters %in% c(4,7,2,17))]


write.table(oligo1,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_1"))
write.table(oligo2,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_2"))

write.table(opc1,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc_1"))
write.table(opc2,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc_2"))

write.table(fibro1,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro_1"))
write.table(fibro2,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro_2"))

