library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
args=commandArgs(trailingOnly=T)
wd=args[1]

#wd="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11cellchat/major"
setwd(wd)
###########obj=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/ONONH_majorclass_40000.rds")
#############obj=subset(obj, subset = majorclass %in% c("Fibroblast","Astrocyte","Oligodendrocyte","Mural_cell","Endothelial_cell","Oligodendrocyte_precursor_cell", "Microglia"))
#data.input <- obj1[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
#labels <- obj1@metadata$majorclass
#meta <- data.frame(labels = labels, row.names = names(labels))

#############Idents(obj) <- obj@meta.data$majorclass

##############obj <- NormalizeData(obj)

##############cellchat <- createCellChat(object = obj, group.by = "ident", assay = "RNA")

###############CellChatDB <- CellChatDB.human

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
################CellChatDB.use <- subsetDB(CellChatDB)
# set the used database in the object
#################cellchat@DB=CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
##############cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 2) # do parallel
options(future.globals.maxSize = 50 * 1024^3)
##############cellchat <- identifyOverExpressedGenes(cellchat)
###########cellchat <- identifyOverExpressedInteractions(cellchat)

########saveRDS(cellchat, file = "cellchat_subset.rds")
cellchat=readRDS("cellchat_subset.rds")
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))
#> [1] 13.20763
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

ptm = Sys.time()
#cellchat <- computeCommunProb(cellchat, type = "triMean")
#> triMean is used for calculating the average gene expression per cell group.
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-02-14 00:32:35.767285]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-02-14 00:33:13.121225]"

#cellchat <- filterCommunication(cellchat, min.cells = 10, min.samples=0)


#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
#df.net <- subsetCommunication(cellchat)

#cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
#cellchat <- aggregateNet(cellchat)
#saveRDS(cellchat, file = "cellchat_subset.rds")

#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))
#> [1] 38.73308

#ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
pdf("major_cellcell_interaction.pdf")
par(mfrow = c(1,2), xpd=TRUE)
p=netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions") #, small.gap = 0.2, big.gap = 2)
p1=netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength") #, small.gap = 0.2, big.gap = 2)
print(p)
print(p1)
dev.off()

mat <- cellchat@net$weight
pdf("major_cellcell_interaction_indi.pdf",width=15,height=15)
par(mfrow = c(5,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]) #, small.gap = 0.2, big.gap = 2)
}
dev.off()

#pdf("major_cellcell_interaction_sign_pair_chord_gene.pdf",width=15,height=15)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
#netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
#netVisual_chord_gene(cellchat,  slot.name = "netP", legend.pos.x = 10)

#dev.off()


pdf("major_cellcell_interaction_sign_pair_bubble.pdf",width=12,height=25)
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
#netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

#######pdf("major_cellcell_interaction_sign_gene_CXCL.pdf",width=12,height=15)
########plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")
#######dev.off()

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
########gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
pdf("major_cellcell_interaction_dominate_sender_receiver_scatter.pdf",width=15)
gg1
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf("major_cellcell_interaction_dominate_sender_receiver_heatmap.pdf",width=15)
ht1 + ht2
dev.off()

library(NMF)
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 2/2
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")

nPatterns = 6
pdf("major_cellcell_interaction_MCP_out_heatmap.pdf",width=10,height=10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()

pdf("major_cellcell_interaction_MCP_out_river.pdf",width=15)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

# dot plot
pdf("major_cellcell_interaction_MCP_out_dot.pdf",width=15)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

selectK(cellchat, pattern = "incoming")

nPatterns = 3

pdf("major_cellcell_interaction_MCP_in_heatmap.pdf", width=10,height=10)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()

# river plot
pdf("major_cellcell_interaction_MCP_in_river.pdf",width=15)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
pdf("major_cellcell_interaction_MCP_in_dot.pdf",width=15)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE, nCores = 1)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf("major_cellcell_interaction_classification_func.pdf",width=15)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE, nCores = 1)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf("major_cellcell_interaction_classification_struc.pdf",width=15)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()
pdf("major_cellcell_interaction_classification_struc_Zoom.pdf",width=15)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 147.8175

saveRDS(cellchat, file = "cellchat_subset.rds")
