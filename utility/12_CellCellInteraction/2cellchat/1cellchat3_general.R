library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
args=commandArgs(trailingOnly=T)
wd=args[1]
#rc_list=read.table(args[2])
sign_pathway=read.table(args[2])
setwd(wd)
future::plan("multisession", workers = 2) # do parallel
options(future.globals.maxSize = 50 * 1024^3)
cellchat=readRDS("cellchat_subset.rds")

# Circle plot
for( path.way in sign_pathway$V1){
pdf(paste0("pathway_",path.way,"_violin.pdf"))
p=plotGeneExpression(cellchat, signaling = path.way , enriched.only = TRUE, type = "violin")
print(p)
dev.off()

pdf(paste0("pathway_",path.way,"_chord_gene.pdf"))
p=netVisual_chord_gene(cellchat, signaling = path.way, legend.pos.x = 8)
print(p)
dev.off()
}


