library(Seurat)
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
args <- commandArgs(trailingOnly = TRUE)
dir=args[1]
name=args[2]
rds=paste0(dir,name,".rds")
meta=paste0(dir,name,".obs.gz")
seurat_obj=readRDS(rds)

sce=as.SingleCellExperiment(seurat_obj)

sce <- slingshot(sce, clusterLabels = 'subclass', reducedDim = 'SCVI')
saveRDS(sce,file=paste0(dir,name,"_slingshot.rds"))
rd <- reducedDims(sce)$SCVI[,1:2]
cl = seurat_obj@meta.data$subclass
lin1 <- getLineages(rd, cl, start.clus = 'OPC_GLIS3+', end.clu="MFOL")
saveRDS(sce,file=paste0(dir,name,"_slingshot.rds"))

rd1=reducedDims(sce)$UMAP[,1:2]
pdf(paste0(dir,name,"_slingshot.pdf"))
palette <- brewer.pal(length(unique(cl)), "Paired")
plot(rd1, col = palette[as.factor(cl)], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black', show.constraints = TRUE)
dev.off()
saveRDS(sce,file=paste0(dir,name,"_slingshot.rds"))
#count_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

#sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

#pto2 <- slingshot(rd2, cl2, omega = TRUE, start.clus = c(1,11))
