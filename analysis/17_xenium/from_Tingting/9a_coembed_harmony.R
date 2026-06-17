# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)

outdir <- "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/harmony"
args <- commandArgs(trailingOnly = TRUE)

#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/neuron_xenium.combined_annot.rds")
#xenium.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/nonneuron_xenium.combined_annot.rds")
#parts <- "nonneuron"

xenium.combined <- args[1]
parts <- args[2]

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

if (parts == "neuron"){
    drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/drg_ref.combined_2k_plus_panel_neuron.rds")
} else if (parts == "nonneuron"){
    drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/drg_ref.combined_2k_plus_panel_nonneuron.rds")
} else if (parts == "all"){
    drg_ref.combined <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_spatial/DRG_Xenium/drg_ref.combined_2k_plus_panel.rds")
}


drg_ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(drg_ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
DefaultAssay(drg_ref.combined)
DefaultAssay(xenium.combined)

drg_ref.combined@assays[["integrated"]]@counts <- drg_ref.combined@assays[["integrated"]]@data
xenium.combined@assays[["integrated"]]@counts <- xenium.combined@assays[["integrated"]]@data

merged_obj <- merge(x = drg_ref.combined, y = xenium.combined)
# merged_obj <- NormalizeData(merged_obj)

features <- rownames(xenium.combined)
features <- intersect(features, rownames(drg_ref.combined))

#merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- Seurat::ScaleData(merged_obj, vars.to.regress = "ident")
merged_obj <- RunPCA(merged_obj, features = features)

library(harmony)

merged_obj <- RunHarmony(merged_obj, "ident")
merged_obj

merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.5, cluster.name = "harmony_clusters")
merged_obj <- FindClusters(merged_obj, resolution = 0.3, cluster.name = "harmony_clusters_p3")
merged_obj <- FindClusters(merged_obj, resolution = 0.6, cluster.name = "harmony_clusters_p8")

merged_obj <- merged_obj %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

pdf(paste0(outdir, "/coembed.combined_harmony_umap_", parts, ".pdf"), width = 10, height = 5)
p1 <- DimPlot(merged_obj, reduction = "umap", group.by = "ident")
p2 <- DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters")
p1 + p2
DimPlot(merged_obj, reduction = "umap", group.by = "Atlas_annotation", label = TRUE,
        repel = TRUE)
DimPlot(merged_obj, reduction = "umap", group.by = "deconv_annot", label = TRUE,
        repel = TRUE)
DimPlot(merged_obj, reduction = "umap", group.by = "final_annot", label = TRUE,
        repel = TRUE)

DimPlot(merged_obj, reduction = "umap", group.by = "spot_class", label = TRUE,
        repel = TRUE)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p3", split.by = "ident", label = T)

#DimPlot(merged_obj, reduction = "umap", group.by = "Atlas_annotation", split.by = "ident", label = T)
dev.off()

table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters)
table(merged_obj$deconv_annot, merged_obj$harmony_clusters)

saveRDS(merged_obj, paste0(outdir, "/coembed.combined_inte_", parts, ".rds"))

#========
merged_obj <- readRDS(paste0(outdir, "/coembed.combined_inte_", parts, ".rds"))
table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters)
table(merged_obj$spot_class, merged_obj$harmony_clusters)
table(merged_obj$deconv_annot, merged_obj$harmony_clusters_p3)
table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters_p3)

DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p3", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p8", label = T)
merged_obj$harmony_clusters_annot <- recode(merged_obj$harmony_clusters, "0" = "Th", "1" = "Calca+Sstr2", "2" = "Trpm8", 
                                            "3" = "Sst", "4" = "Mrgprd", "5" = "Mrgpra3+Trpv1", "6" = "Calca+Smr2", 
                                            "7" = "Calca+Bmpr1b", "8" = "Calca+Oprk1", "9" = "Mrgprd", "10" = "Ntrk3high+Ntrk2", "11" = "?", 
                                            "12" = "Pvalb", "13" = "..", "14" = "Calca+Dcn", "15" = "Pvalb", "16" = "Ntrk3low+Ntrk2",
                                            "17" = "Calca+Adra2a", "18" = "Rxfp1", "19" = "Mrgpra3+Mrgprb4", "20" = "...", 
                                            "21" = "Ntrk3high+S100a16", "22" = "Ntrk3high+Ntrk2", "23" = "Calca+Oprk1", "24" = "Pvalb", "25" = "Mrgpra3+Trpv1", 
                                            "26" = "Atf3", "27" = "Mrgpra3+Trpv1")
table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters_annot)

doublet_df <- merged_obj@meta.data[merged_obj@meta.data$ident == "Xenium",]
doublet_df <- doublet_df[doublet_df$spot_class == "doublet_certain",]
table(doublet_df$deconv_annot, doublet_df$harmony_clusters)


