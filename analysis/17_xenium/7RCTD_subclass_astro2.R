# Usage:
# $ Rscript /dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_atlas/TG_atlas/3_nonneuron_RCTD.R"sampleid"

# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)

args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(R_SCRIPT = "/pub/junw42/local/miniconda/envs/spatial1/bin/Rscript")
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/spacexr/")
dir.create(output_PATH, recursive = TRUE)



xenium_data_list=readRDS(paste0(output_PATH1, target_tissue,"_8_xenium_data_list_comb_cutoff20.rds"))
label=args[2]

#tg.ref=readRDS(args[1])
clu=args[6]
i=as.numeric(args[7])
ref.combined=readRDS(args[1]) #readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table(args[4], header=T,sep=",",row.names=1) #read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)

ref.combined <- AddMetaData(ref.combined, metadata = meta)

#######ref.combined = subset(ref.combined, subset = sampleid == args[5])

tg.ref=ref.combined

features <- rownames(xenium_data_list[[1]])
features <- intersect(features, rownames(tg.ref))


#table(tg.ref$subclass)
#table(tg.ref$subclass1)

#for(i in 3:length(xenium_data_list)){
#for(i in 3:5){
t=unique(xenium_data_list[[i]]@meta.data$slide_id)
bc=read.table(file=paste0(args[3],"_",t,"_",label),header=T,row.names=1)

#bc=read.table(file=paste0(args[3],"_",t,"_",label),header=T,row.names=1, sep=",")
xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)

#xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$barcode)
xenium_data_list[[i]]=xenium_data_list[[i]][features,]
xenium_obj=xenium_data_list[[i]]


 
  
  query.counts <- GetAssayData(xenium_obj, assay = "Xenium", slot = "counts")
  coords <- xenium_obj@images$fov$centroids@coords
  #coords <- xenium_obj@meta.data[, c("x", "y")]
  rownames(coords) <- rownames(xenium_obj@meta.data)
  coords <- as.data.frame(coords)
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
  
#Idents(tg.ref) <- "subclass"
#Idents(tg.ref) <- "subclass1"

#  Idents(tg.ref) <- "subclass2"
  Idents(tg.ref) <- clu
#  Idents(tg.ref) <- "subclass"

  counts <- GetAssayData(tg.ref, assay = "RNA", slot = "counts")
#  cluster <- as.factor(tg.ref[,clu])
##cluster <- as.factor(tg.ref$subclass)

#cluster <- as.factor(tg.ref$subclass1)

  cluster <- as.factor(tg.ref@meta.data[[clu]])


  #  cluster <- as.factor(tg.ref$subclass2)
#    cluster <- as.factor(tg.ref$subclass)

  names(cluster) <- colnames(tg.ref)
  nUMI <- colSums(counts)
  names(nUMI) <- colnames(tg.ref)
  reference <- Reference(counts, cluster, nUMI)
  
  # run RCTD with many cores
#  RCTD <- create.RCTD(query, reference, max_cores = 28, UMI_min = 0, gene_cutoff = 0.000125, fc_cutoff = 0.25, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.5) #nonneuron1
#  RCTD <- create.RCTD(query, reference, max_cores = 28, UMI_min = 0, gene_cutoff = 0.0001, fc_cutoff = 0.15, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.3) #nonneuron1

  #  RCTD <- create.RCTD(query, reference, max_cores = 28, UMI_min = 0, gene_cutoff = 0.000125, fc_cutoff = 0.25, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.5) #nonneuron
 #RCTD <- create.RCTD(query, reference, max_cores = 28, UMI_min = 0)
# RCTD <- create.RCTD(query, reference, max_cores = 12, UMI_min = 0)
 RCTD <- create.RCTD(query, reference, max_cores = 1, UMI_min = 0, gene_cutoff = 0.000125, fc_cutoff = 0.20, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.3)

Sys.setenv(PATH = paste("/pub/junw42/local/miniconda/envs/spatial1/bin/", Sys.getenv("PATH"), sep=":"))
stopifnot(nzchar(Sys.which("Rscript")))
system("Rscript --version")
#RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
Sys.which("R")
Sys.which("Rscript")
#class(param)
#param

  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
#  saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_subclass2_RCTD_results_nonneuron.rds"))
#  saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_", clu, "_RCTD_results_nonneuron.rds"))
  saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_", clu, "_RCTD_results_nonneuron020.rds"))

######  saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_RCTD_results_nonneuron.rds"))
#   saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron.rds"))
#   saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron025.rds"))
#   saveRDS(RCTD, paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron015.rds"))

#}
