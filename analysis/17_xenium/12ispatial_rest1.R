library(future)
plan("multicore")  # or "multisession" if you're on Windows
options(future.globals.maxSize = 80 * 1024^3)  # 50 GB

library(Seurat)
library(ggplot2)
library(dplyr)
library(FNN)
library(harmony)
source("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/from_Tingting/ispatial/iSpatial.R")


args <- commandArgs(trailingOnly = TRUE)

# read xenium data
#slide_name <- "Round2_Slide04"
slide_name = "ONH_1" #<- args[1]
#print(slide_name)
outdir <- paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/imputation/",slide_name, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]

wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1", ".rds"))

genes_select = rownames(merged_obj)

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

integrated=readRDS(paste0(outdir, "/", slide_name, "_integrated.rds"))

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
DimPlot(integrated_scRNA, group.by = "majorclass", label = T, repel = T)
DimPlot(integrated_merFISH, group.by = "harmony_anno", label = T, repel = T)
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
  }, mc.cores = 8)
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

saveRDS(enhancer_expr, paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))


enhancer_expr = do.call(cbind, enhancer_expr)
colnames(enhancer_expr) = colnames(integrated_merFISH)

saveRDS(enhancer_expr, paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))


#enhancer_expr <- readRDS(paste0(outdir, "/", slide_name, "_enhancer_expr.rds"))


# assign inferred expression value to new assay
storage.mode(enhancer_expr) <- "numeric"
rownames(enhancer_expr)=rownames(integrated_merFISH@assays$RNA@data)
integrated_merFISH[[infered.assay]] = SeuratObject::CreateAssayObject(data = as(enhancer_expr, "dgCMatrix"))

DefaultAssay(integrated_merFISH) <- infered.assay

integrated_merFISH = Seurat::DietSeurat(integrated_merFISH, assays = infered.assay)

saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH.rds"))

#integrated_merFISH <- readRDS(paste0(outdir, "/", slide_name, "_integrated_merFISH.rds"))

# add spatial information
integrated_merFISH@images = spRNA_images


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

############

library(pheatmap)
library(RColorBrewer)
library(viridis)

get_gene_xenium <- function(xenium.combined, gene, type){
  expr_mat <- xenium.combined@assays$enhanced$data[rownames(xenium.combined@assays$enhanced$data) %in% gene, xenium.combined@meta.data$harmony_anno == type]
  return(mean(expr_mat))
}

marker_list=c("GAD1","SLC4A11","GRM6","ARR3","PTPRB","BICC1","ONECUT1","CD69","F13A1", "MLANA","RGR","CD74","NOTCH3", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ")

#marker_list <- c("Atf3", "Calca", "Adra2a", "Bmpr1b",
#                 "Oprk1", "Smr2", "Sstr2", "Mrgpra3", "Mrgprb4", "Trpv1", "Mrgprd",
#                 "Ntrk3", "S100a16", "Ntrk2", "Sst","Th",
#                 "Trpm8")
#types <- sort(unique(integrated_merFISH@meta.data$harmony_anno))
types=c("AC","Astrocyte","BC","Cone","Endothelial-cell","Fibroblast","HC","Immune-cell","Macrophage","Melanocyte","MG","Microglia","Mural-cell","Oligodendrocyte","Oligodendrocyte-precursor-cell","RGC","RPE","Rod", "Schwann-cell")
mean_list <- list()
for (i in c(1:length(marker_list))){
  gene <- marker_list[i]
  mean_list[[i]] <- sapply(types, function(x){get_gene_xenium(integrated_merFISH, gene, x)})
}
hm_df <- as.data.frame(do.call(cbind, mean_list))
colnames(hm_df) <- marker_list
rownames(hm_df) <- types
saveRDS(hm_df, paste0(outdir, "/", slide_name, "_hm_df.rds"))

pdf(paste0(outdir, "/", slide_name, "_imputed_marker_heatmap.pdf"), width = 8, height = 6)
p <- pheatmap(hm_df, scale = "column", cluster_cols = F, cluster_rows = F, color = viridis(n = 10, alpha = 1,
                                                                                           begin = 0, end = 1, option = "viridis"), border_color = "black")
print(p)
dev.off()

#############
integrated_merFISH <- FindVariableFeatures(integrated_merFISH, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(integrated_merFISH)
integrated_merFISH <- ScaleData(integrated_merFISH, features = all.genes)
integrated_merFISH <- RunPCA(integrated_merFISH, npcs = 30, features = VariableFeatures(integrated_merFISH))
integrated_merFISH <- RunUMAP(integrated_merFISH, reduction = "pca", dims = 1:30)
integrated_merFISH <- FindNeighbors(integrated_merFISH, dims = 1:30, reduction = "pca")
integrated_merFISH <- FindClusters(integrated_merFISH, resolution = 1, cluster.name = "imputed_w1_by_slide")

saveRDS(integrated_merFISH, paste0(outdir, "/", slide_name, "_integrated_merFISH_clustering.rds"))

pdf(paste0(outdir, "/", slide_name, "_umap.pdf"), width = 8, height = 5)
DimPlot(integrated_merFISH)
DimPlot(integrated_merFISH, group.by = "harmony_anno", label = T, repel = T)
dev.off()
#pdf(paste0(outdir, "/", slide_name, "_neuron_umap_tglabel.pdf"), width = 10, height = 5)
#DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T, split.by = "tglabel")
#dev.off()
#pdf(paste0(outdir, "/", slide_name, "_neuron_umap_branchlabel.pdf"), width = 15, height = 5)
#DimPlot(integrated_merFISH, group.by = "final_annot", label = T, repel = T, split.by = "branchlabel")
#dev.off()



