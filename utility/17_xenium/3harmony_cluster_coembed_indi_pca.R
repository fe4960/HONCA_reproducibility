# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)

####outdir <- "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/harmony"
args <- commandArgs(trailingOnly = TRUE)

#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/neuron_xenium.combined_annot.rds")
#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/nonneuron_xenium.combined_annot.rds")
#parts <- "nonneuron"


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

#xenium.combined <- args[1]
#parts <- args[2]

#============================================================================================================
# co-embedding - Harmony

#xenium.combined <- readRDS("DRG_Xenium/xenium.combined.rds")
#drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_spatial/DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")

#non_neuron_types <- c("Endothelial", "Fibroblast", "Immune", "Pericyte", "Satellite glia", "Schwann_M", "Schwann_N")
#drg_ref.combined@meta.data$is_neuron <- ifelse(drg_ref.combined@meta.data$Atlas_annotation %in% non_neuron_types, 0, 1)

#drg_ref.combined <- subset(drg_ref.combined, subset = is_neuron == 1)
#drg_ref.combined <- subset(drg_ref.combined, subset = is_neuron == 0)
#dim(drg_ref.combined)

#saveRDS(drg_ref.combined, "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/drg_ref.combined_2k_plus_panel_neuron.rds")
#saveRDS(drg_ref.combined, "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/drg_ref.combined_2k_plus_panel_nonneuron.rds")

#if (parts == "neuron"){
#    drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/drg_ref.combined_2k_plus_panel_neuron.rds")
#} else if (parts == "nonneuron"){
#    drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/drg_ref.combined_2k_plus_panel_nonneuron.rds")
#} else if (parts == "all"){
#    drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_spatial/DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")
#}

#ref.combined=readRDS(paste0(output_PATH, target_tissue,"_ref.combined_5k_plus_panel.rds"))

ref.combined=readRDS(args[1]) #readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table(args[2], header=T,sep=",",row.names=1) #read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)

ref.combined <- AddMetaData(ref.combined, metadata = meta)

ref.combined = subset(ref.combined, subset = sampleid == args[5])



#####process xenium data

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_cutoff20.rds"))
label=args[4]

features <- rownames(xenium_data_list[[1]])
features <- intersect(features, rownames(ref.combined))


for(i in length(xenium_data_list)){

bc=read.table(file=paste0(args[3],"_",i),header=T,row.names=1)
xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)
xenium_data_list[[i]]=xenium_data_list[[i]][features,]
}





xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")
xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)

DefaultAssay(xenium.combined) <- "integrated"
######
ref.combined = ref.combined[features,]
ref.combined <- NormalizeData(ref.combined)
ref.combined <- ScaleData(ref.combined, features = features, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, features = features, verbose = FALSE)

obj_list <- list("ref" = ref.combined, "xenium" = xenium.combined)

#features <- rownames(xenium.combined)
#features <- intersect(features, rownames(ref.combined))

coembed.anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca")
coembed.combined <- IntegrateData(anchorset = coembed.anchors, new.assay.name = "coembed")
DefaultAssay(coembed.combined) <- "coembed"



#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(ref.combined)
DefaultAssay(xenium.combined)

######for multiple sample
#ref.combined@assays[["integrated"]]@counts <- ref.combined@assays[["integrated"]]@data
#####for one sample

#new_assay <- CreateAssayObject(counts = ref.combined@assays[["RNA"]]@data)

#ref.combined@assays[["integrated"]]=new_assay
#Key(ref.combined[["integrated"]])="integrated_"
#ref.combined@assays[["integrated"]]@counts <- ref.combined@assays[["RNA"]]@data

#DefaultAssay(ref.combined)="integrated"

xenium.combined@assays[["integrated"]]@counts <- xenium.combined@assays[["integrated"]]@data

merged_obj <- merge(x = ref.combined, y = xenium.combined)
# merged_obj <- NormalizeData(merged_obj)

features <- rownames(xenium.combined)
features <- intersect(features, rownames(ref.combined))

#merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- Seurat::ScaleData(merged_obj, vars.to.regress = "ident")
merged_obj <- RunPCA(merged_obj, features = features)

library(harmony)

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

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_sampleid_",label,".pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="sampleid")+ NoLegend()
dev.off()

merged_obj@meta.data$slide="snRNA"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_1"]="xenium_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_2"]="xenium_2"

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

#####

