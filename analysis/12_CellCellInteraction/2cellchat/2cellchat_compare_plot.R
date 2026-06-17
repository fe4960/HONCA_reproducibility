load(paste0(wd2, "cellchat_merged_drg_g",g,"_chip_sham.RData")

pdf(paste0(wd2,"compareInt_circle_macro.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)


netVisual_diffInteraction(cellchat, weight.scale = T, sources.use ="Macrophage")
netVisual_diffInteraction(cellchat, weight.scale = T, targets.use ="Macrophage")
dev.off()


pdf(paste0(wd2,"compareInt_circle_macro_strength.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)


netVisual_diffInteraction(cellchat, weight.scale = T, sources.use ="Macrophage", measure = "weight", arrow.size=1)
netVisual_diffInteraction(cellchat, weight.scale = T, targets.use ="Macrophage", measure = "weight", arrow.size=1)
dev.off()

pdf(paste0(wd2,"netVisual_bubble_macro_source.pdf"),width=12,height=12)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use = "Macrophage", targets.use = NULL ,  comparison = c(1, 2), angle.x = 45, max.dataset = 2)

netVisual_bubble(cellchat, sources.use = "Macrophage", targets.use = NULL ,  comparison = c(1, 2), angle.x = 45, max.dataset = 1)
dev.off()


pdf(paste0(wd2,"netVisual_bubble_macro_target.pdf"),width=12,height=12)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_bubble(cellchat, sources.use = NULL, targets.use = "Macrophage" ,  comparison = c(1, 2), angle.x = 45, max.dataset = 2)

netVisual_bubble(cellchat, sources.use = NULL, targets.use = "Macrophage" ,  comparison = c(1, 2), angle.x = 45, max.dataset = 1)
dev.off()


#netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Macrophage", targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Macrophage", targets.use = NULL, stacked = F, do.stat = TRUE)
pdf(paste0(wd2,"compareInt_circle_macro_rankNet.pdf"),width=12,height=6)
gg1 + gg2
dev.off()


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL , targets.use = "Macrophage", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL , targets.use = "Macrophage", stacked = F, do.stat = TRUE)
pdf(paste0(wd2,"compareInt_circle_macro_rankNet_target.pdf"),width=12,height=6)
gg1 + gg2
dev.off()


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Macrophage", targets.use = "Neuron", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Neuron", targets.use = "Macrophage", stacked = F, do.stat = TRUE)
pdf(paste0(wd2,"compareInt_circle_macro_neuron.pdf"),width=12,height=6)
 gg2
dev.off()

library(ComplexHeatmap)

load(paste0(wd2, "cellchat_object.drg_g",g,"_list.RData"))

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
pdf(paste0(wd2,"compareInt_signalingRole_heatmap.pdf"),width=12,height=6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


##########

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "chip"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
####net.up <- subsetCommunication(cellchat, net = net, datasets = "chip",ligand.logFC = 0.05, receptor.logFC = NULL)
net.up <- subsetCommunication(cellchat, net = net, datasets = "chip",ligand.logFC = 0.02, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
####net.down <- subsetCommunication(cellchat, net = net, datasets = "sham",ligand.logFC = -0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "sham",ligand.logFC = -0.02, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

#df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")

pdf(paste0(wd2,"compareInt_heatmap_macro.pdf"),width=12,height=6)
# Chord diagram
par(mfrow = c(2,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = "Macrophage", targets.use = NULL, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[2]], sources.use = "NULL", targets.use = "Macrophage", slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = "Neuron", targets.use = "Macrophage", slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

#netVisual_chord_gene(object.list[[1]], sources.use = NULL, targets.use = "Macrophage", slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway
dev.off()


##############
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = "Macrophage", targets.use = c("Neuron","SatGlia","Schwann_M"), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = "Neuron", targets.use = c("Neuron","SatGlia","Schwann_M"), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pdf(paste0(wd2,"netVisual_bubble_neuron.pdf"),width=12,height=6)
gg1 + gg2
dev.off()


###########

#Identify signaling groups based on their functional similarity
#ptm = Sys.time()

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf(paste0(wd2,"netVisual_pairwise.pdf"),width=12,height=6)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()
#> 2D visualization of signaling networks from datasets 1 2

###############


pathways.show <- c("ADGRL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste0(wd2,"netVisual_aggregate_ADGRL.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()


cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("sham", "chip")) # set factor level

#cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("sham", "chip")) # set factor level
pdf(paste0(wd2,"netVisual_geneExp.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
plotGeneExpression(cellchat, signaling = "ADGRL", split.by = "datasets", colors.ggplot = T, type = "violin")
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")
dev.off()

#############
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Neuron", signaling.exclude = NULL)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophage", signaling.exclude = NULL)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
pdf(paste0(wd2,"netAnalysis_signalingChanges.pdf"),width=12,height=6)
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


#############
library(ggplot2)
pdf(paste0(wd2,"netVisual_chord_gene_macro_staglia.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = "Macrophage", targets.use = "Schwann_M", slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))+scale_color_brewer(palette="Set3")
netVisual_chord_gene(object.list[[2]], sources.use = "Schwann_M", targets.use = "Macrophage", slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))+scale_color_brewer(palette="Set3")
dev.off()
#netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
########
pdf(paste0(wd2,"netVisual_chord_gene_macro_staglia_002.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = "Macrophage", targets.use = "SatGlia", slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))+scale_color_brewer(palette="Set3")
#####netVisual_chord_gene(object.list[[2]], sources.use = "SatGlia", targets.use = "Macrophage", slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))+scale_color_brewer(palette="Set3")
dev.off()

##########

pathways.show <- c("CXCL")
pdf(paste0(wd2,"netVisual_chord_cell_CXCL.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:1) { #1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,  title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
dev.off()


pathways.show <- c("APP")
pdf(paste0(wd2,"netVisual_chord_cell_APP.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,  title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
dev.off()
