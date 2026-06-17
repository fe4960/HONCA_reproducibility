library(Seurat)
library(patchwork)

slide=c("ONH_1", "ONH_2", "ONH_3", "PP_1", "PP_2")
#gene=""
#BMP7
for( t in slide){
fn=paste0("HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/imputation/", t, "/", t, "_integrated_merFISH.rds")
obj=readRDS(fn)
#cell=rownames(obj@meta.data[obj$harmony_anno=="Astrocyte",])

cell <- WhichCells(obj, expression = harmony_anno == "Astrocyte")

p1 <- ImageFeaturePlot(obj, features = c("APOE","SOX9","SOX2"), cells =cell)
p2 <- ImageDimPlot(obj, molecules =c("APOE","SOX9","SOX2") , nmols = 10000, alpha = 0.3, mols.cols = "red", group.by="harmony_anno", cells=cell)
pdf(paste0("HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/imputation/", t, "/", t, "_integrated_astro.pdf"), width = 20, height = 15)
#pdf(paste0("HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/imputation/", t, "/", t, "_integrated_BMP7.pdf"), width = 18, height = 9)
print(p1 + p2)
dev.off()
}
