cellchat=readRDS("cellchat_subset.rds")
pathways.show <- c("BMP")
pdf("major_cellcell_interaction_BMP_circle.pdf")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()


pdf("major_cellcell_interaction_BMP_chord.pdf")
par(mfrow=c(1,1))
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()


pdf("major_cellcell_interaction_BMP_heatmap.pdf")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

pdf("major_cellcell_interaction_BMP_contr.pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
