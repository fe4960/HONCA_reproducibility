# slurmucihpc.sh -m 50G -p standard -a ruic20_lab -- Rscript test_ispatial_p3.R slide_name
library(future)
plan("multicore")  # or "multisession" if you're on Windows
options(future.globals.maxSize = 200 * 1024^3)  # 100 GB

# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(FNN)

args <- commandArgs(trailingOnly = TRUE)

# read xenium data
#slide_name <- "Round2_Slide03"
slide_name <- args[1]
print(slide_name)
outdir <- paste0("/dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_atlas/analysis/imputation/", slide_name, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

integrated <- readRDS(paste0(outdir, "/", slide_name, "_integrated.rds"))
enhancer_expr <- readRDS(paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))


integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30, reduction = "harmony")


DefaultAssay(integrated) <- "RNA"
#integrated = suppressWarnings(Seurat::DietSeurat(integrated, counts = F, assays = "RNA"))
integrated_scRNA = suppressWarnings(subset(integrated, subset = tech == "scRNA"))
integrated_merFISH = suppressWarnings(subset(integrated, subset = tech == "scRNA", invert = TRUE))
cells_name = colnames(integrated_scRNA)

pdf(paste0(outdir, "/", slide_name, "_ispatial_inte.pdf"), width = 8, height = 5)
DimPlot(integrated, group.by = "tech")
DimPlot(integrated_scRNA, group.by = "Atlas_annotation", label = T, repel = T)
DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T)
dev.off()

#===


infered.assay <- "enhanced"

# assign inferred expression value to new assay
integrated_merFISH[[infered.assay]] = SeuratObject::CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))

DefaultAssay(integrated_merFISH) <- infered.assay

integrated_merFISH = Seurat::DietSeurat(integrated_merFISH, assays = infered.assay)

saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH.rds"))



#=====

integrated_merFISH <- FindVariableFeatures(integrated_merFISH, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated_merFISH)
integrated_merFISH <- ScaleData(integrated_merFISH, features = all.genes)
integrated_merFISH <- RunPCA(integrated_merFISH, npcs = 30, features = VariableFeatures(integrated_merFISH))
integrated_merFISH <- RunUMAP(integrated_merFISH, reduction = "pca", dims = 1:30)
integrated_merFISH <- FindNeighbors(integrated_merFISH, dims = 1:30, reduction = "pca")
integrated_merFISH <- FindClusters(integrated_merFISH, resolution = 0.3, cluster.name = "imputed_p3_by_slide")

saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH_clustering.rds"))

pdf(paste0(outdir, "/", slide_name, "_neuron_umap.pdf"), width = 8, height = 5)
DimPlot(integrated_merFISH)
DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T)
dev.off()
pdf(paste0(outdir, "/", slide_name, "_neuron_umap_tglabel.pdf"), width = 10, height = 5)
DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T, split.by = "tglabel")
dev.off()
pdf(paste0(outdir, "/", slide_name, "_neuron_umap_branchlabel.pdf"), width = 15, height = 5)
DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T, split.by = "branchlabel")
dev.off()

#====

