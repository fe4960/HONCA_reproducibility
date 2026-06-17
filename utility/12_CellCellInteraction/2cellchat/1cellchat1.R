library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
wd="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11cellchat/major"
setwd(wd)
obj=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/ONONH_majorclass_40000.rds")
obj=subset(obj, subset = majorclass %in% c("Fibroblast","Astrocyte","Oligodendrocyte","Mural_cell","Endothelial_cell","Oligodendrocyte_precursor_cell", "Microglia"))
#data.input <- obj1[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
#labels <- obj1@metadata$majorclass
#meta <- data.frame(labels = labels, row.names = names(labels))

Idents(obj) <- obj@meta.data$majorclass

obj <- NormalizeData(obj)

cellchat <- createCellChat(object = obj, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.human

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
# set the used database in the object
cellchat@DB=CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 2) # do parallel
options(future.globals.maxSize = 50 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

saveRDS(cellchat, file = "cellchat_subset.rds")

