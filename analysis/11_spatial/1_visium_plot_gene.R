library(Seurat)
library(ggplot2)
library(dplyr)
#source("/dfs3b/ruic20_lab/singlecell/tingtiny/fetal_visiumHD/functions_visiumHD.R")
source("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/10_spatial/0functions_visiumHD.R")
# pre-processed 10w fetal VisiumHD data
object <- readRDS("/dfs3b/ruic20_lab/singlecell/tingtiny/fetal_visiumHD/D1_10w/analysis/object.rds")
# pre-processed 15w fetal VisiumHD data
object <- readRDS("/dfs3b/ruic20_lab/singlecell/tingtiny/fetal_visiumHD/A1_15w/object_deconvolution.rds")

main="/dfs3b/ruic20_lab/junw42/"

time="10w"

cell="Fibro"

gene="SMOC2"
# 1) to plot the gene expression on the spatial
#p=plot_gene_VHD(object, "BEST1")
p=plot_gene_VHD(object, gene)

pdf(paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_plot.pdf"), width=6, height=7)
print(p)
dev.off()
# 2) to highlight the spots with the gene expressed on the spatial
#plot_gene_VHD_highlight(object, "BEST1") # highlight the bins 
plot_gene_VHD_highlight(object, gene) # highlight the bins 

pdf(paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_highlight.pdf"), width=6, height=7)
print(p)
dev.off()


# 3) To highlight the spots with the gene1 expressed, gene2 expressed, and both genes co-expressed on the spatial separately
#plot_gene_VHD_highlight_2genes(object, "BEST1", "RPE")
gene1="SHOX"
p=plot_gene_VHD_highlight_2genes(object, gene, gene1)

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_",gene1,"_highlight.png"), plot = p, width = 11, height = 5, dpi = 300)

#pdf(paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_",gene1,"_highlight.pdf"), width=15, height=7)
#print(p)
#dev.off()
gene="MGP"
gene1="APOD"
p=plot_gene_VHD_highlight_2genes(object, gene, gene1)

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_",gene1,"_highlight.png"), plot = p, width = 11, height = 5, dpi = 300)

gene="NOX4"

p=plot_gene_VHD(object, gene)

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_sclera_highlight.png"), plot = p, width = 6, height = 5, dpi = 300)

gene="CEMIP"
gene1="PKP2"
p=plot_gene_VHD_highlight_2genes(object, gene, gene1)

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_",gene1,"_arachnoid_highlight.png"), plot = p, width = 11, height = 5, dpi = 300)

gene="SLC47A1"
gene1="THSD4"

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_",gene1,"_dura_highlight.png"), plot = p, width = 11, height = 5, dpi = 300)


gene="AL445250.1"
gene1="ABCA9"

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_",gene1,"_pia_highlight.png"), plot = p, width = 11, height = 5, dpi = 300)


gene="BMP5"

p=plot_gene_VHD(object, gene)

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_rpe_highlight.png"), plot = p, width = 6, height = 5, dpi = 300)


gene="SLC22A3"

p=plot_gene_VHD(object, gene)

ggsave(filename = paste0(main,"HCA_ON/data/8_spatial/visiumHD/",time,"_",cell,"_",gene,"_x_highlight.png"), plot = p, width = 6, height = 5, dpi = 300)

