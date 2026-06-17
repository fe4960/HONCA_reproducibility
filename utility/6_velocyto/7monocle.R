library(Seurat)
library(monocle3)
args <- commandArgs(trailingOnly = TRUE)
dir=args[1]
name=args[2]
rds=paste0(dir,name,".rds")
meta=paste0(dir,name,".obs.gz")
seurat_obj=readRDS(rds)
count_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cell_metadata=read.table(meta,sep=",",header=T, row.names=1)
gene=data.frame(gene_short_name=rownames(count_matrix))
rownames(gene)=rownames(count_matrix)
cds <- new_cell_data_set(count_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "sampleid")
#, residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "subclass")


cds <- cluster_cells(cds)


plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "subclass",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#####

root_cells <- colnames(cds)[cds$subclass == "OPC_GLIS3+"]
cds <- order_cells(cds, root_cells = root_cells)

## Step 7: Save rds
saveRDS(cds, paste0(dir, name, "_monocel3_cds.rds"))

write.csv(pseudotime(cds, reduction_method = "UMAP"), paste0(dir,name, "_cds_pseudotime_monocle3.csv"))

## Step 8: QC
png(paste0(dir,"/figures/", name, "_cds_pseudotime_monocle3_subclass.png"), width = 1000, height = 1000)
plot_cells(cds, color_cells_by = "subclass", cell_size = 1)
dev.off()


png(paste0(dir,"/figures/", name, "_cds_pseudotime_monocle3_pseudotime.pdf"), width = 1000, height = 1000)
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1)
dev.off()


mk <- c("SOX10","OLIG1","OLIG2","PLP1","MYRF","MOG","MAG","PTPRZ1","PDGFRA","VCAN","GPR17","BCAS1","PTGDS")

png(paste0(dir,"/figures/", name, "_cds_pseudotime_monocle3_subclass_mk.png"), width = 1000, height = 1000)
plot_cells(cds,genes=mk, label_cell_groups=FALSE,show_trajectory_graph=FALSE)
dev.off()

df <- data.frame(cds$subclass, pseudotime(cds, reduction_method = "UMAP"))
colnames(df) <- c("celltype", "time")
#df$days <- factor(df$days)
#df$region <- plyr::mapvalues(x = df$region, from = c("FR", "NR"), to = c("Macula",   "Peripheral"))
p <- ggplot(df, aes(x = celltype, y = time, fill = celltype)) + geom_boxplot(outlier.shape = NA) +
    xlab("Cell type") + ylab("Pseudotime") + scale_fill_brewer(palette="Paired") + theme(text = element_text(size = 25))
pdf(paste0(dir,"/figures/", nam, "_cds_pseudotime_monocle3_pseudotime_subclass.png"), width = 10, height = 5)
print(p)
dev.off()
