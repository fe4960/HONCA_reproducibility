# Usage:
# $ Rscript /dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_atlas/TG_atlas/3_nonneuron_RCTD.R"sampleid"

# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)

args <- commandArgs(trailingOnly = TRUE)

# read xenium data
#sampleid <- "Slide178_Region_1"
sampleid <- args[1]
outdir <- paste0("/dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_atlas/analysis/TG_Xenium_Section_process/", sampleid, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(paste0(outdir, sampleid, "_RCTD_results_nonneuron.rds"))) {
  
  xenium_obj <- readRDS(paste0(outdir, sampleid, "_nonneuron_obj.rds"))
  
  
  query.counts <- GetAssayData(xenium_obj, assay = "Xenium", slot = "counts")
  coords <- xenium_obj@images$fov$centroids@coords
  #coords <- xenium_obj@meta.data[, c("x", "y")]
  rownames(coords) <- rownames(xenium_obj@meta.data)
  coords <- as.data.frame(coords)
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
  
  
  tg.ref <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_ref/TG_snRNA_yang_2022_downsample5k.rds")
  table(tg.ref@meta.data$Atlas_annotation)
  neuron_types <- c("Mrgprd", "Ntrk3high+Ntrk2", "Trpm8",           
                    "Calca+Sstr2", "Th", "Ntrk3low+Ntrk2",     
                    "Ntrk3high+S100a16", "Mrgpra3+Mrgprb4", "Calca+Smr2",       
                    "Calca+Bmpr1b", "Sst","Calca+Oprk1",        
                    "Calca+Adra2a", "Mrgpra3+Trpv1", "Atf3")               
  
  nonneuron_types <- c("Oligodendrocyte-like", "Arachnoid", "Endothelial",   
                       "Immune", "Dura", "Satellite glia",   
                       "Schwann_N", "Pia", "Pericyte",         
                       "Astrocyte-like", "Schwann_M")
  
  tg.ref$majorclass <- ifelse(tg.ref$Atlas_annotation %in% neuron_types, "Neuron", tg.ref@meta.data$Atlas_annotation)
  table(tg.ref$majorclass)
  
  #tg.ref <- subset(tg.ref, subset = Atlas_annotation != "Atf3")
  #table(tg.ref@meta.data$Atlas_annotation)
  Idents(tg.ref) <- "majorclass"
  
  counts <- GetAssayData(tg.ref, assay = "RNA", slot = "counts")
  cluster <- as.factor(tg.ref$majorclass)
  names(cluster) <- colnames(tg.ref)
  nUMI <- colSums(counts)
  names(nUMI) <- colnames(tg.ref)
  reference <- Reference(counts, cluster, nUMI)
  
  # run RCTD with many cores
  RCTD <- create.RCTD(query, reference, max_cores = 28, UMI_min = 0)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  saveRDS(RCTD, paste0(outdir, sampleid, "_RCTD_results_nonneuron.rds"))
  
}
