# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)


od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="5" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue1, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
#path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0)

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

ref.combined=readRDS(paste0(output_PATH1, target_tissue1,"_ref.combined_5k_plus_panel.rds"))
xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(ref.combined)
DefaultAssay(xenium.combined)

ref.combined@assays[["integrated"]]@counts <- ref.combined@assays[["integrated"]]@data
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

pdf(paste0(output_PATH, target_tissue, "coembed.combined_harmony_umap", ".pdf"), width = 12, height = 5)
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

saveRDS(merged_obj, paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", ".rds"))

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", ".rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct.pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="majorclass",label=T)+ NoLegend()
dev.off()


merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_1"] <- "PP_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_2"] <- "PP_2"

merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_1"] <- "ONH_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_2"] <- "ONH_2"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_3"] <- "ONH_3"



pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_sampleid.pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="sampleid")+ NoLegend()
dev.off()

merged_obj@meta.data$slide="snRNA"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_1"]="xenium_PP_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_2"]="xenium_PP_2"

merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_1"]="xenium_ONH_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_2"]="xenium_ONH_2"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_3"]="xenium_ONH_3"


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_slide.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", split.by="slide")+ NoLegend()
dev.off()

table(merged_obj@meta.data[,c("majorclass","harmony_clusters_p6")])

tb=table(merged_obj@meta.data[,c("harmony_clusters_p6","majorclass")])
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_heatmap.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(p6_clu=rownames(tb),ct=max_colnames)

merged_obj@meta.data$harmony_anno="Unk"

for( i in anno$p6_clu){
merged_obj@meta.data$harmony_anno[merged_obj@meta.data$harmony_clusters_p6==i]=anno[anno$p6_clu==i,]$ct
}

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_slide_ct.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno", split.by="slide")+ NoLegend()
dev.off()

saveRDS(merged_obj, paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony", ".rds"))

merged_obj@meta.data[merged_obj@meta.data$harmony_anno=="T_cell",]$harmony_anno="Immune_cell"
merged_obj@meta.data[merged_obj@meta.data$harmony_anno=="Oligodendrocyte_precursor_cell",]$harmony_anno="OPC"

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ident_ct_anno.pdf"), width = 10, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno", split.by="ident", label =TRUE) + NoLegend()
dev.off()


#####

#========
######merged_obj <- readRDS(paste0(outdir, "/coembed.combined_inte_", parts, ".rds"))
#table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters)
#table(merged_obj$spot_class, merged_obj$harmony_clusters)
#table(merged_obj$deconv_annot, merged_obj$harmony_clusters_p3)
#table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters_p3)

################
######DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p3", label = T)
#######DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters", label = T)
#######DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p8", label = T)
########merged_obj$harmony_clusters_annot <- recode(merged_obj$harmony_clusters, "0" = "Th", "1" = "Calca+Sstr2", "2" = "Trpm8", 
#####                                            "3" = "Sst", "4" = "Mrgprd", "5" = "Mrgpra3+Trpv1", "6" = "Calca+Smr2", 
####                                            "7" = "Calca+Bmpr1b", "8" = "Calca+Oprk1", "9" = "Mrgprd", "10" = "Ntrk3high+Ntrk2", "11" = "?", 
#####                                            "12" = "Pvalb", "13" = "..", "14" = "Calca+Dcn", "15" = "Pvalb", "16" = "Ntrk3low+Ntrk2",
#####                                            "17" = "Calca+Adra2a", "18" = "Rxfp1", "19" = "Mrgpra3+Mrgprb4", "20" = "...", 
####                                            "21" = "Ntrk3high+S100a16", "22" = "Ntrk3high+Ntrk2", "23" = "Calca+Oprk1", "24" = "Pvalb", "25" = "Mrgpra3+Trpv1", 
#####                                            "26" = "Atf3", "27" = "Mrgpra3+Trpv1")
#####table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters_annot)

#####doublet_df <- merged_obj@meta.data[merged_obj@meta.data$ident == "Xenium",]
#####doublet_df <- doublet_df[doublet_df$spot_class == "doublet_certain",]
#####table(doublet_df$deconv_annot, doublet_df$harmony_clusters)


