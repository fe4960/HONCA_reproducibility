library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
args=commandArgs(trailingOnly=T)
wd=args[1]
wd1=args[2]
wd2=paste0(wd1,"disease_normal_cellchat_comparison/")
#wd2=paste0(wd1,"chip_sham_cellchat_comparison/")
dir.create(wd2,recursive = TRUE)
g=args[3]
#wd="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11cellchat/major"
#setwd(wd2)
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
#cellchat=readRDS("cellchat_subset.rds")

#ptm = Sys.time()

#####cellchat.sham <- readRDS(paste0(wd,"cellchat_subset.rds"))
#######cellchat.chip <- readRDS(paste0(wd1,"cellchat_subset.rds"))

#cellchat.sham <- readRDS(paste0(wd,g, "_cellchat_subset.rds"))
#cellchat.chip <- readRDS(paste0(wd1,g, "_cellchat_subset.rds"))



#######object.list <- list(sham = cellchat.sham, chip = cellchat.chip)
########cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))


# Users can now export the merged CellChat object and the list of the two separate objects for later use
#save(object.list, file = paste0(wd2, "cellchat_object.g_",g,"_list.RData"))
#save(cellchat, file = paste0(wd2, "cellchat_merged_g_",g,"_chip_sham.RData"))

#######save(object.list, file = paste0(wd2,"ellchat_object.g_",g,"_list.RData"))
######save(cellchat, file = paste0(wd2, "cellchat_merged_g_",g,"_chip_sham.RData"))

load(paste0(wd2,"cellchat_object.g_",g,"_list.RData"))
load(paste0(wd2, "cellchat_merged_g_",g,"_chip_sham.RData"))


####gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
######gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
######pdf(paste0(wd2,"g_", g, "_compareInt_bar.pdf"),width=12,height=6)
######gg1 + gg2
#####dev.off()

#####pdf(paste0(wd2,"g_", g, "_compareInt_circle.pdf"),width=12,height=6)
######par(mfrow = c(1,2), xpd=TRUE)
#####netVisual_diffInteraction(cellchat, weight.scale = T)
######netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#######3dev.off()

########gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
#######gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
#######pdf(paste0(wd2,"g_", g, "_compareInt_heatmap.pdf"),width=12,height=6)
########gg1 + gg2
########dev.off()


########
#############num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
##############weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
###########gg <- list()
##########for (i in 1:length(object.list)) {
#########  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
###########}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

############pdf(paste0(wd2,"g_", g, "_compareInt_signalRole_scatter.pdf"),width=12,height=6)
##########patchwork::wrap_plots(plots = gg)
#########dev.off()
#########
####gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OLIGO2_LRRC7+")

#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophage", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
#####gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OLIGO1")

#gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Neuron", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
#gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Schwann_M", signaling.exclude = c("MIF"))
#######gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Oligodendrocyte_precursor_cell")
#gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OLIGO1_SVEP1+")
#gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OLIGO1_RBFOX1+")



#########3pdf(paste0(wd2,"g_", g, "_compareInt_signalChanges_scatter_Oligo.pdf"),width=25,height=10)
#########patchwork::wrap_plots(plots = list(gg1,gg2,gg3))

#patchwork::wrap_plots(plots = list(gg1,gg2,gg3, gg4, gg5))
##########dev.off()


#Astro_ONHON1_SLC4A11+MARCH1+
#OLIGO1_SVEP1+
#Fibroblast
#Astro_ONH1_GABBR2+SV2B+
#Astro_ONH2_PAX5+GABBR2+
#Microglia
#Endothelial_cell
#Oligodendrocyte_precursor_cell
#OLIGO2_LRRC7+
#Mural_cell
#OLIGO1
#Astro_ONHON2_CST3+APOE+
#Astro_retina1_PAX5+ME1+
#OLIGO1_RBFOX1+

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ON")+theme(
    axis.title = element_text(size = 25),
    axis.text  = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25),
    plot.title = element_text(size = 25)
  )


#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ON1_WNK2+TSHZ2+")
#gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ON2_ACTN1+SERPINA3+")
#gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ON3_NR4A3+FOS+")

gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ONHON")+theme(
    axis.title = element_text(size = 25),
    axis.text  = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25),
    plot.title = element_text(size = 25)
  )

#####gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ONHON2_CST3+APOE+")

#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophage", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
#gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ONH1_GABBR2+SV2B+")
gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ONH")+theme(
    axis.title = element_text(size = 25),
    axis.text  = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25),
    plot.title = element_text(size = 25)
  )


#gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Neuron", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
#gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Schwann_M", signaling.exclude = c("MIF"))
#gg7 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_ONH2_PAX5+GABBR2+")
#gg8 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_retina1_PAX5+ME1+")
gg8 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astro_retina")+theme(
    axis.title = element_text(size = 25),
    axis.text  = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25),
    plot.title = element_text(size = 25)
  )



pdf(paste0(wd2,"g_", g, "_compareInt_signalChanges_scatter_astro.pdf"),width=25,height=10)
patchwork::wrap_plots(plots = list(gg1, gg4,  gg6,  gg8))

#patchwork::wrap_plots(plots = list(gg1,gg2,gg3, gg4, gg5, gg6, gg7, gg8))
dev.off()
