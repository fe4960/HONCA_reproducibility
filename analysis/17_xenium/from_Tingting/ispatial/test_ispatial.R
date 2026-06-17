# slurmucihpc.sh -m 150G -p standard -a ruic20_lab -- Rscript test_ispatial.R slide_name
library(future)
plan("multicore")  # or "multisession" if you're on Windows
options(future.globals.maxSize = 200 * 1024^3)  # 100 GB

# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(FNN)

source("/dfs3b/ruic20_lab/tingty7/codes/Xenium/iSpatial.R")

args <- commandArgs(trailingOnly = TRUE)

# read xenium data
#slide_name <- "Round2_Slide04"
slide_name <- args[1]
print(slide_name)
outdir <- paste0("/dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_atlas/analysis/imputation/", slide_name, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


dir_PATH <- "/dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_atlas/analysis/refine_annot/TG_Xenium_Section_process"
Round_Slide_pattern <- paste0("^", slide_name, ".*_xenium_obj_cutoff20_final_annot_niche\\.rds$")

obj_PATHs <- list.files(path = dir_PATH, pattern = Round_Slide_pattern, recursive = TRUE, full.names = TRUE)

xenium_data_list <- list()
for (i in c(1:length(obj_PATHs))){
  obj_PATH <- obj_PATHs[i]
  print(obj_PATH)
  xenium_data <- readRDS(obj_PATH)
  xenium_data <- subset(xenium_data, subset = majorclass == "Neuron")
  xenium_data <- subset(xenium_data, subset = spot_class == "singlet")
  xenium_data_list[[i]] <- xenium_data
}

print("merging data...")
xenium.merge <- merge(x = xenium_data_list[[1]], y = xenium_data_list[2:length(xenium_data_list)])
xenium.merge <- JoinLayers(xenium.merge)

xenium.merge <- NormalizeData(xenium.merge)
features <- rownames(xenium.merge)
xenium.merge <- ScaleData(xenium.merge, features = features)
xenium.merge <- RunPCA(xenium.merge, npcs = 30, features = features)
xenium.merge <- RunUMAP(xenium.merge, reduction = "pca", dims = 1:30)
xenium.merge <- FindNeighbors(xenium.merge, dims = 1:30, reduction = "pca")
xenium.merge <- FindClusters(xenium.merge, resolution = 0.3, cluster.name = "xenium_clusters_p3_by_slide")

dim(xenium.merge) # 

pdf(paste0(outdir, "/", slide_name, "_neuron_xenium_umap_annot.pdf"), width = 7, height = 5)
DimPlot(xenium.merge, label = T)
DimPlot(xenium.merge, label = T, group.by = "final_annot", repel = T)
dev.off()

saveRDS(xenium.merge, paste0(outdir, "/", slide_name, "_neuron_xenium.merge.rds"))

xenium.combined <- xenium.merge
#xenium.combined <- readRDS(paste0(outdir, "/", slide_name, "_neuron_xenium.merge.rds"))
rm(xenium.merge)

options(Seurat.object.assay.version = "v3")
spRNA <- CreateSeuratObject(xenium.combined@assays$Xenium$counts, assay = "RNA", meta.data = xenium.combined@meta.data)
spRNA <- NormalizeData(spRNA)

tg.ref_full <- readRDS("/dfs3b/ruic20_lab/tingty7/projects/TG_spatial/TG_ref/TG_neuron_snRNA_yang_2022.rds")
table(tg.ref_full@meta.data$Atlas_annotation)

options(Seurat.object.assay.version = "v3")
scRNA <- CreateSeuratObject(tg.ref_full@assays$RNA$counts, assay = "RNA", meta.data = tg.ref_full@meta.data)
scRNA <- NormalizeData(scRNA)

# parameters

dims = 1:30
k.neighbor = 30
infered.assay = "enhanced"
weighted.KNN = TRUE
RNA.weight = 0.5
n.core = 8
correct.spRNA = FALSE
correct.scRNA = FALSE
correct.weight.NN = 0.2
correct.neighbor = 5

spRNA$tech = "spatial"
scRNA$tech = "scRNA"

if(is.null(spRNA@assays$RNA)){
  stop(paste(spRNA, " do not have 'RNA' assay."))
}

if(is.null(scRNA@assays$RNA)){
  stop(paste(scRNA, " do not have 'RNA' assay."))
}

if(length(spRNA@assays$RNA@data) == 0){
  stop(paste(spRNA, " is not normlized. Run Seurat::NormalizeData."))
}

if(length(scRNA@assays$RNA@data) == 0){
  stop(paste(scRNA, " is not normlized. Run Seurat::NormalizeData."))
}

if(correct.spRNA){
  message("Stablize spatial transcriptome.")
  spRNA = stabilize_expr(spRNA, 
                         neighbor = correct.neighbor,
                         weight.NN = correct.weight.NN,
                         n.core = n.core, npcs = length(dims))
}

if(correct.scRNA){
  message("Stablize single cell RNAseq.")
  scRNA = stabilize_expr(scRNA, 
                         neighbor = correct.neighbor, 
                         weight.NN = correct.weight.NN,
                         n.core = n.core, npcs = length(dims))
}

genes_select = intersect(rownames(scRNA), rownames(spRNA))

if(length(genes_select) < 10){
  stop("Too few intesected genes between scRNA and spRNA.")
}

# run pca
spRNA = Seurat::FindVariableFeatures(spRNA, verbose = FALSE)
SeuratObject::VariableFeatures(spRNA) = SeuratObject::VariableFeatures(spRNA)[SeuratObject::VariableFeatures(spRNA) %in% genes_select]
spRNA = Seurat::ScaleData(spRNA, verbose = FALSE)
spRNA = Seurat::RunPCA(spRNA, npcs = length(dims), verbose = FALSE)

SeuratObject::VariableFeatures(scRNA) = SeuratObject::VariableFeatures(spRNA)
scRNA = Seurat::ScaleData(scRNA, verbose = FALSE)
scRNA = Seurat::RunPCA(scRNA, npcs = length(dims), verbose = FALSE)


# merge two objects
message("1st level integration")
anchors <- Seurat::FindIntegrationAnchors(object.list = list(spRNA, scRNA), 
                                          normalization.method = "LogNormalize", 
                                          reduction = "rpca", 
                                          anchor.features = SeuratObject::VariableFeatures(spRNA),
                                          k.anchor = 20,
                                          dims = dims,  
                                          verbose = FALSE)
integrated <- Seurat::IntegrateData(anchorset = anchors, 
                                    dims = dims, 
                                    normalization.method = "LogNormalize",
                                    verbose = TRUE)
rm(anchors)
#gc()

# count level normalization
message("normalization")
norm_data = expm1(integrated@assays$RNA@data)


# trim high expression gene, which affect the mean of genes
trim_quantil1 = Matrix::rowMeans(scRNA@assays$RNA@data[genes_select, ])
trim_quantil1 = trim_quantil1[trim_quantil1 < quantile(trim_quantil1, prob = 0.98) ] 

trim_quantil2 = Matrix::rowMeans(spRNA@assays$RNA@data[genes_select, ])
trim_quantil2 = trim_quantil2[trim_quantil2 < quantile(trim_quantil2, prob = 0.98) ] 

trim_genes = unique(c(names(trim_quantil1), names(trim_quantil2)))

# normalize factor based on trimmed genes
norm_factor = Matrix::colMeans(norm_data[trim_genes, ]) 
norm_factor = as.numeric(norm_factor)
norm_factor = norm_factor/mean(norm_factor)

norm_data@x <- norm_data@x / rep.int(norm_factor, diff(norm_data@p))
integrated@assays$RNA@data = log1p(norm_data)

# remove cell with number of expressed gene < 98%
integrated = integrated[, norm_factor != 0]

#spRNA_images = suppressWarnings(spRNA@images)
spRNA_images = spRNA@meta.data[, c("x", "y")]

#rm(norm_data, norm_factor, scRNA, spRNA)
rm(norm_data, norm_factor)
#gc()


# check normalization
#avg_expr = Seurat::AverageExpression(integrated, group.by="tech", slot="data")
#boxplot(log1p(avg_expr$RNA[genes_select, ]))

# 2nd level integration
message("2nd level integration")
#SeuratObject::VariableFeatures(integrated, assay = "integrated") = genes_select
integrated = Seurat::ScaleData(integrated, vars.to.regress = "tech", #features = genes_select,
                               assay = "integrated", 
                               verbose = FALSE)
integrated = Seurat::RunPCA(integrated, 
                            npcs = length(dims), 
                            assay = "integrated", 
                            verbose = FALSE)

integrated <- suppressWarnings( harmony::RunHarmony(
  object = integrated,
  group.by.vars = 'tech',
  plot_convergence = F,
  theta = 3,
  lambda = 0.5,
  assay.use = "integrated",
  verbose = FALSE,
  max.iter.harmony = 20
))

saveRDS(integrated, paste0(outdir, "/", slide_name, "_integrated.rds"))

#integrated <- readRDS(paste0(outdir, "/", slide_name, "_integrated.rds"))
# ignore...==========
# find neighbors
#integrated = Seurat::FindNeighbors(integrated, 
#                                   k.param = k.neighbor,
#                                   reduction="harmony",
#                                   dims=dims, 
#                                   return.neighbor = T, 
#                                   assay = "integrated")

#neigbors = suppressWarnings(sapply(colnames(subset(integrated, subset = tech == "scRNA", invert = TRUE)), function(i){
#  Seurat::TopNeighbors(integrated@neighbors$integrated.nn, i, n = k.neighbor)
#}))

#neigbors = as.data.frame(neigbors)

#DefaultAssay(integrated) <- "RNA"
#integrated = suppressWarnings(Seurat::DietSeurat(integrated, counts = F, assays = "RNA"))
#integrated_scRNA = suppressWarnings(subset(integrated, subset = tech == "scRNA"))
#integrated_merFISH = suppressWarnings(subset(integrated, subset = tech == "scRNA", invert = TRUE))
#cells_name = colnames(integrated_scRNA)

#neigbors = as.list(neigbors)

# only select neighbor in scRNA data
#neigbors = lapply(neigbors, function(x){
#  x = x[x %in% cells_name]
#})


##### modify by TY#####

## --- settings from your script ---
k.use   <- k.neighbor
dims.use <- dims  # e.g., 1:30

## --- split cells by tech ---
is_sc <- integrated$tech == "scRNA"
is_sp <- !is_sc

## --- Harmony embeddings (shared co-embedding) ---
H     <- Embeddings(integrated, reduction = "harmony")[, dims.use, drop = FALSE]
H_sc  <- H[is_sc, , drop = FALSE]
H_sp  <- H[is_sp, , drop = FALSE]

sc_names <- rownames(H_sc)
sp_names <- rownames(H_sp)

## guard: k cannot exceed # scRNA cells
k.use <- min(k.use, nrow(H_sc))

## --- KNN: for each spatial cell, find k nearest scRNA neighbors ---
kn <- FNN::get.knnx(data = H_sc, query = H_sp, k = k.use)

## --- SAME OUTPUT STRUCTURE as your original 'neigbors' (list) ---
neigbors <- lapply(seq_len(nrow(H_sp)), function(i) sc_names[kn$nn.index[i, ]])
names(neigbors) <- sp_names

#====

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
enhancer_expr = integrated_merFISH@assays$RNA@data

saveRDS(enhancer_expr, paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))

# infer expression via scRNA
message("infer expression.")
if(weighted.KNN){
  # dynamics proportion to assign scRNA values to merFISH
  enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
    cell_neighbors = neigbors[[cell]]
    if (length(cell_neighbors) == 0){
      integrated@assays$RNA@data[, cell]
    }else{
      cor_dist = sparse.cor(integrated@assays$RNA@data[genes_select, c(cell, cell_neighbors)])[,1]
      cor_dist[is.na(cor_dist)] <- 0 
      cor_dist[cor_dist < 0] <- 0
      cor_dist = cor_dist ** 2
      # normalized correlation distance matrix
      cor_dist = cor_dist / sum(cor_dist)
      cor_dist = c((1-RNA.weight) * cor_dist[1], RNA.weight * cor_dist[-1]) # cor_dist[1] is the cell in spRNA, here = 1
      
      # inner produce
      infer_expr = integrated@assays$RNA@data[, c(cell, cell_neighbors)] %*% cor_dist 
      infer_expr[,1] 
    }
  }, mc.cores = n.core)
}else{
  # constant proportion to assign scRNA values to merFISH
  enhancer_expr = parallel::mclapply(colnames(enhancer_expr), function(cell){
    cell_neighbors = neigbors[[cell]]
    if (length(cell_neighbors) == 0){
      enhancer_expr[, cell]
    }else if (length(cell_neighbors) == 1){
      (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * integrated_scRNA@assays$RNA@data[,cell_neighbors]
    }else{
      (1-RNA.weight) * enhancer_expr[, cell] + RNA.weight * Matrix::rowMeans(integrated_scRNA@assays$RNA@data[,cell_neighbors])
    }
  }, mc.cores = n.core)
}

enhancer_expr = do.call(cbind, enhancer_expr)
colnames(enhancer_expr) = colnames(integrated_merFISH)

saveRDS(enhancer_expr, paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))

#enhancer_expr <- readRDS(paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))


# assign inferred expression value to new assay
integrated_merFISH[[infered.assay]] = SeuratObject::CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))

DefaultAssay(integrated_merFISH) <- infered.assay

integrated_merFISH = Seurat::DietSeurat(integrated_merFISH, assays = infered.assay)

saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH.rds"))

#integrated_merFISH <- readRDS(paste0(outdir, "/", slide_name, "_integrated_merFISH.rds"))

# add spatial information
integrated_merFISH@images = spRNA_images

#for(img in names(integrated_merFISH@images)){
#  integrated_merFISH@images[[img]]@coordinates <- integrated_merFISH@images[[img]]@coordinates[
#    rownames(integrated_merFISH@images[[img]]@coordinates) %in% colnames(integrated_merFISH), ]
#}

integrated_merFISH@images <- integrated_merFISH@images[rownames(integrated_merFISH@images) %in% colnames(integrated_merFISH), ]

saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH.rds"))

#========
dim(integrated_merFISH)


cor_list <- as.numeric()
i <- 1
for (gene in genes_select){
 sp_exp <- spRNA@assays$RNA$counts[gene,] 
 infer_exp <- integrated_merFISH@assays$enhanced$data[gene, names(sp_exp)]
 cor_list[[i]] <- cor(sp_exp, infer_exp)
 i <- i + 1
}

summary(cor_list)
saveRDS(cor_list, paste0(outdir, "/", slide_name, "cor_list.rds"))
#cor_list <- readRDS(paste0(outdir, "/", slide_name, "cor_list.rds"))


# stats
m  <- mean(cor_list, na.rm = TRUE)
md <- median(cor_list, na.rm = TRUE)

# for placing labels slightly above the curve
d     <- density(cor_list, na.rm = TRUE)
y_top <- max(d$y, na.rm = TRUE)

p <- ggplot(data.frame(x = cor_list), aes(cor_list)) +
  geom_density(fill = "grey90", color = "grey35", linewidth = 1) +
  geom_rug(alpha = 0.2) +
  # dashed lines
  geom_vline(xintercept = m,  linetype = "dashed", linewidth = 0.7, color = "#1f77b4") +
  geom_vline(xintercept = md, linetype = "dashed", linewidth = 0.7, color = "#d62728") +
  # text labels
  annotate("text", x = m,  y = y_top * 1.05, label = sprintf("Mean = %.3f", m),
           color = "#1f77b4", hjust = 0.5, vjust = 0) +
  annotate("text", x = md, y = y_top * 1.15, label = sprintf("Median = %.3f", md),
           color = "#d62728", hjust = 0.5, vjust = 0) +
  labs(x = "PCC", y = "Density", title = paste0("PCC - ", slide_name)) +
  coord_cartesian(ylim = c(0, y_top * 1.25)) +
  theme_classic()

pdf(paste0(outdir, "/", slide_name, "_pcc_distribution.pdf"), width = 4, height = 3)
print(p)
dev.off()

#=========

library(pheatmap)
library(RColorBrewer)
library(viridis)

get_gene_xenium <- function(xenium.combined, gene, type){
  expr_mat <- xenium.combined@assays$enhanced$data[rownames(xenium.combined@assays$enhanced$data) %in% gene, xenium.combined@meta.data$final_annot == type]
  return(mean(expr_mat))
}

#get_gene_xenium(integrated_merFISH, "Dcn", "Atf3")
#marker_list <- c("Fabp7", "Pecam1", "Dcn", "Mpz", "Ptprc", "Notch3", "Scn7a")
marker_list <- c("Atf3", "Calca", "Adra2a", "Bmpr1b", 
                 "Oprk1", "Smr2", "Sstr2", "Mrgpra3", "Mrgprb4", "Trpv1", "Mrgprd",
                 "Ntrk3", "S100a16", "Ntrk2", "Sst","Th",
                 "Trpm8")
types <- sort(unique(integrated_merFISH@meta.data$final_annot))
mean_list <- list()
for (i in c(1:length(marker_list))){
  gene <- marker_list[i]
  mean_list[[i]] <- sapply(types, function(x){get_gene_xenium(integrated_merFISH, gene, x)})
}
hm_df <- as.data.frame(do.call(cbind, mean_list))
colnames(hm_df) <- marker_list
rownames(hm_df) <- types
#hm_df <- hm_df[types,]
#hm_df <- hm_df[, c("Atf3", "Calca", "Adra2a", "Bmpr1b", "Dcn",
#                   "Oprk1", "Smr2", "Sstr2", "Mrgpra3", "Mrgprb4", "Trpv1", "Mrgprd",
#                   "Ntrk3", "S100a16", "Ntrk2","Pvalb", "Rxfp1","Sst","Th",
#                   "Trpm8")]

#hm_df[hm_df > 2.5] <- 2.5
#hm_df[hm_df < -2.5] <- -2.5
saveRDS(hm_df, paste0(outdir, "/", slide_name, "_neuron_hm_df.rds"))



#=======
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


###==
pdf(paste0(outdir, "/", slide_name, "_neuron_imputed_marker_heatmap.pdf"), width = 8, height = 6)
p <- pheatmap(hm_df, scale = "column", cluster_cols = F, cluster_rows = F, color = viridis(n = 10, alpha = 1, 
                                                                                           begin = 0, end = 1, option = "viridis"), border_color = "black")
print(p)
dev.off()
