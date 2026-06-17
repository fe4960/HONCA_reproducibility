# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)
library(harmony)

options(future.globals.maxSize = 50 * 1024^3)

####outdir <- "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/harmony"
args <- commandArgs(trailingOnly = TRUE)

#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/neuron_xenium.combined_annot.rds")
#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/nonneuron_xenium.combined_annot.rds")
#parts <- "nonneuron"


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
# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0)

#xenium.combined <- args[1]
#parts <- args[2]

#============================================================================================================
# co-embedding - Harmony


ref.combined=readRDS(args[1]) #readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")


#####process xenium data

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))
label=args[2]

features <- rownames(xenium_data_list[[1]])
features <- intersect(features, rownames(ref.combined))


for(i in 1:length(xenium_data_list)){
t=unique(xenium_data_list[[i]]@meta.data$slide_id)
bc=read.table(file=paste0(args[3],"_",t,"_",label),header=T,row.names=1)
xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)
xenium_data_list[[i]]=xenium_data_list[[i]][features,]
}





xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "cca")

#xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")
xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)

DefaultAssay(xenium.combined) <- "integrated"
######

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(ref.combined)="integrated"
DefaultAssay(xenium.combined)="integrated"

######for multiple sample
#ref.combined@assays[["integrated"]]@counts <- ref.combined@assays[["integrated"]]@data
#####for one sample

#######new_assay <- CreateAssayObject(counts = ref.combined@assays[["RNA"]]@data)

#########ref.combined@assays[["integrated"]]=new_assay
#########Key(ref.combined[["integrated"]])="integrated_"
ref.combined@assays[["integrated"]]@counts <- ref.combined@assays[["integrated"]]@data

########DefaultAssay(ref.combined)="integrated"

xenium.combined@assays[["integrated"]]@counts <- xenium.combined@assays[["integrated"]]@data

merged_obj <- merge(x = ref.combined, y = xenium.combined)
# merged_obj <- NormalizeData(merged_obj)

features <- rownames(xenium.combined)
features <- intersect(features, rownames(ref.combined))

#merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- Seurat::ScaleData(merged_obj, vars.to.regress = "ident")
merged_obj <- RunPCA(merged_obj, features = features)



merged_obj <- RunHarmony(merged_obj, "ident")
merged_obj

merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.5, cluster.name = "harmony_clusters_p5")
merged_obj <- FindClusters(merged_obj, resolution = 0.3, cluster.name = "harmony_clusters_p3")
merged_obj <- FindClusters(merged_obj, resolution = 0.6, cluster.name = "harmony_clusters_p6")

merged_obj <- merged_obj %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

pdf(paste0(output_PATH, target_tissue, "coembed.combined_harmony_umap_",label, ".pdf"), width = 12, height = 5)
####p1 <- DimPlot(merged_obj, reduction = "umap", group.by = "ident")
####p2 <- DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters")
####p1 + p2
###DimPlot(merged_obj, reduction = "umap", group.by = "Atlas_annotation", label = TRUE,repel = TRUE)
####DimPlot(merged_obj, reduction = "umap", group.by = "deconv_annot", label = TRUE, repel = TRUE)
####DimPlot(merged_obj, reduction = "umap", group.by = "final_annot", label = TRUE,repel = TRUE)

######DimPlot(merged_obj, reduction = "umap", group.by = "spot_class", label = TRUE,repel = TRUE)
#DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p3", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p5", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p6", split.by = "ident", label = T)

#DimPlot(merged_obj, reduction = "umap", group.by = "Atlas_annotation", split.by = "ident", label = T)
dev.off()

#table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters)
#table(merged_obj$deconv_annot, merged_obj$harmony_clusters)

saveRDS(merged_obj, paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, ".rds"))

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, ".rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,".pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="subclass",label=T,split.by="ident")+ NoLegend()
dev.off()


merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_1"] <- "PP_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_2"] <- "PP_2"

merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_1"] <- "ONH_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_2"] <- "ONH_2"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_3"] <- "ONH_3"



pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_sampleid_",label,".pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="sampleid")+ NoLegend()
dev.off()

merged_obj@meta.data$slide="snRNA"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_1"]="xenium_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_2"]="xenium_2"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_1"] <- "xenium_3"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_2"] <- "xenium_4"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_3"] <- "xenium_5"


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_slide_",label,".pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", split.by="slide")+ NoLegend()
dev.off()

table(merged_obj@meta.data[,c("majorclass","harmony_clusters_p3")])

tb=table(merged_obj@meta.data[,c("harmony_clusters_p3","subclass")])
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_heatmap_",label,".pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(p3_clu=rownames(tb),ct=max_colnames)

merged_obj@meta.data$harmony_anno="Unk"

for( i in anno$p3_clu){
merged_obj@meta.data$harmony_anno[merged_obj@meta.data$harmony_clusters_p3==i]=anno[anno$p3_clu==i,]$ct
}
saveRDS(merged_obj, paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, ".rds"))

#####
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_anno.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno",label=F,split.by="slide") #+ NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_subclass.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="subclass",label=F,split.by="slide") #+ NoLegend()
dev.off()


