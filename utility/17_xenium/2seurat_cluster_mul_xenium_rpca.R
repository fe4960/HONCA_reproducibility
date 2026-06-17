library(Seurat)
library(ggplot2)
library(dplyr)
#library(arrow)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="5" #args[6]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
#dir.create(output_PATH, recursive = TRUE)
# read xenium data

cut=as.numeric(args[1])

#xenium_data_list_comb=readRDS(paste0(output_PATH, target_tissue, "_xenium_data_list_comb_cutoff20.rds"))
########
#####process xenium data

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_6_xenium_data_list_comb_cutoff20.rds"))
label=args[2]
fc=as.numeric(args[4])
features <- rownames(xenium_data_list[[1]])
#features <- intersect(features, rownames(ref.combined))

#quantile(merged_obj@meta.data[merged_obj@meta.data$slide!="snRNA","nCount_Xenium"])
#0%2125%5050%8275%121100%487
#quantile(merged_obj@meta.data[merged_obj@meta.data$slide!="snRNA","nFeature_Xenium"])
#0%325%2550%3475%44100%105
for(i in 1:length(xenium_data_list)){
t=unique(xenium_data_list[[i]]@meta.data$slide_id)
bc=read.table(file=paste0(args[3],"_",t,"_",label),header=T,row.names=1)
#xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = (barcode %in% bc$x) & (nCount_Xenium >=cut) & (nFeature_Xenium >=fc) )

xenium_data_list[[i]] <- xenium_data_list[[i]][,
  ((xenium_data_list[[i]]$barcode %in% bc$x) &
  (xenium_data_list[[i]]$nCount_Xenium >= cut) &
  (xenium_data_list[[i]]$nFeature_Xenium >= fc)) 
]


#xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = barcode %in% bc$x)
#xenium_data_list[[i]]=xenium_data_list[[i]][features,]
colnames(xenium_data_list[[i]])=paste0(t,"_",colnames(xenium_data_list[[i]]))

}


xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")

#xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "cca")
xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)

DefaultAssay(xenium.combined) <- "integrated"


########
#xenium_data_list_comb <- lapply(X = xenium_data_list_comb, FUN = function(x) {
#  x <- NormalizeData(x)
#})

#features <- rownames(xenium_data_list_comb[[1]])
#xenium_data_list_comb <- lapply(X = xenium_data_list_comb, FUN = function(x) {
#  x <- ScaleData(x, features = features, verbose = FALSE)
#  x <- RunPCA(x, features = features, verbose = FALSE)
#})

#xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list_comb, anchor.features = features, reduction = "rpca")
#xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)

#DefaultAssay(xenium.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
xenium.combined <- ScaleData(xenium.combined, verbose = FALSE)
xenium.combined <- RunPCA(xenium.combined, npcs = 30, verbose = FALSE)
xenium.combined <- RunUMAP(xenium.combined, reduction = "pca", dims = 1:30)
xenium.combined <- FindNeighbors(xenium.combined, reduction = "pca", dims = 1:30)
xenium.combined <- FindClusters(xenium.combined, resolution = 0.3, cluster.name = "xenium_clusters_p3")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.5, cluster.name = "xenium_clusters_p5")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.8, cluster.name = "xenium_clusters_p8")
xenium.combined <- FindClusters(xenium.combined, resolution = 0.2, cluster.name = "xenium_clusters_p2")
xenium.combined <- FindClusters(xenium.combined, resolution = 1.5, cluster.name = "xenium_clusters_15")
xenium.combined <- FindClusters(xenium.combined, resolution = 1, cluster.name = "xenium_clusters_1")

dim(xenium.combined) # >20: TM 818; CB 36134; surface: 28505 ONONH: 478 49290

#####xenium.combined=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium.combined_cutoff50.rds")
#xenium.combined <- FindClusters(xenium.combined, resolution = 0.1, cluster.name = "xenium_clusters_p1")
saveRDS(xenium.combined, paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))



#####pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_umap_1_cutoff",cut,"_fc_",fc,"_p1_6.pdf"), width = 6, height = 5)
#######DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p1") + ggtitle("Xenium Clusters 1")
#DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p5") + ggtitle("Xenium Clusters p5")

#######dev.off()

##########xenium.combined@meta.data$xenium_clusters_p2_anno=xenium.combined@meta.data$xenium_clusters_p2

#########xenium.combined@meta.data[xenium.combined@meta.data$xenium_clusters_p2==4,]$xenium_clusters_p2_anno=3
########xenium.combined@meta.data[xenium.combined@meta.data$xenium_clusters_p2==5,]$xenium_clusters_p2_anno=3
#######xenium.combined@meta.data[xenium.combined@meta.data$xenium_clusters_p2==6,]$xenium_clusters_p2_anno=0

#saveRDS(xenium.combined, paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium.combined_cutoff50.rds"))

##############pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_umap_1_cutoff50_p2.pdf"), width = 6, height = 5)
############DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p2_anno") + ggtitle("Xenium Clusters p2 anno")+scale_color_brewer(palette="Set3")
#DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p5") + ggtitle("Xenium Clusters p5")

###########dev.off()

pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_umap_p2_cutoff",cut,"_fc_",fc,"_6.pdf"), width = 6, height = 5)
#DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p1" ) + ggtitle("Xenium Clusters p1")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p2" ) + ggtitle("Xenium Clusters p2")

DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p3") + ggtitle("Xenium Clusters p3")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p5") + ggtitle("Xenium Clusters p5")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p8") + ggtitle("Xenium Clusters p8")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_15") + ggtitle("Xenium Clusters 1.5")
DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_1") + ggtitle("Xenium Clusters 1")

DimPlot(xenium.combined, group.by = "slide_id") + ggtitle("Sections")
FeaturePlot(xenium.combined, label = T, features = "nCount_Xenium", max.cutoff = "q90")
FeaturePlot(xenium.combined, label = T, features = "nFeature_Xenium")
FeaturePlot(xenium.combined, label = T, features = "cell_size")
dev.off()

#pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_umap_1_cutoff",cut,"_fc_",fc,".pdf"), width = 6, height = 5)
#dev.off()

markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue,"_",label,"_cutoff",cut,"_fc_",fc, "_markers_1_6.csv"))

#table(xenium.combined@meta.data$xenium_clusters_p3)

#markers <- FindAllMarkers(xenium.combined, only.pos = TRUE)
#markers <- markers %>%
#  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 1)

#write.csv(markers, paste0(output_PATH, target_tissue,"_",label,"_cutoff",cut,"_fc_",fc, "_markers_p3.csv"))


#markers <- FindAllMarkers(xenium.combined, group.by = "xenium_clusters_p5", only.pos = TRUE)
#markers <- markers %>%
#  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 1)

#write.csv(markers, paste0(output_PATH, target_tissue,"_",label,"_cutoff",cut,"_fc_",fc, "_markers_p5.csv"))


#marker_genes <- c("MLANA", "SCN7A", "PECAM1", "CCL21", "PTPRC", "ACTA2", "DES", "MME", "POU6F2", "KRT5", "PAX6", "RGS5", "PDGFRB")
################marker_genes=c("GFAP","AC092957.1","MECOM","VWF","BICC1","GABRG3","HKDC1","F13A1","PAX3","KCNQ3","ADAM28","MYH11","CARMN","RNF220","CTNNA3","NXPH1","CSMD1","BCL11B","PYHIN1","CD247","KIT","SKAP1","BANK1")
###################pdf(paste0(output_PATH, target_tissue, "_xenium_vln_dot_p5_cutoff",cut,".pdf"), width = 12, height = 5)
################VlnPlot(xenium.combined, features = marker_genes, pt.size = 0)
###############DotPlot(xenium.combined, features = marker_genes)
#############dev.off()
#saveRDS(xenium.combined, paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,".rds"))


#markers=read.table( paste0(output_PATH, target_tissue, "_markers_p2.csv"))

#########marker_genes=c("THSD7B","TENM1","SPRR1B","CYSRT1","ATP8A2","VIM","AQP4","FAM78B","NELL2","EFEMP1","SCGB2A1","CD74","PTPRC","C1QC","TLR2","CEMIP","LUM","IL13RA2","CLMP", "MMP2")

#########pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_dot_p2_cutoff",cut,".pdf"), width = 12, height = 5)
########DotPlot(xenium.combined, features = marker_genes)
########dev.off()



#opc=c("CSMD1","TNR","OLIG1","OLIG2","OPCML","SNTG1")

#peri=c("RGS5","ABCC9","TRPC4")

#oligo=c("MOG","MOBP","ST18")

#mk=c(opc,peri, oligo)

#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_vln_dot_p2_cutoff",cut,"_fc_",fc,"_6.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = mk, pt.size = 0)
DotPlot(xenium.combined, features = mk)
dev.off()


pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_vln_cellsize_p2_cutoff",cut,"_fc_",fc,"_6.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = "cell_size", pt.size = 0)
#DotPlot(xenium.combined, features = mk)
dev.off()


pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_vln_nCount_Xenium_p2_cutoff",cut,"_fc_",fc,"_6.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = "nCount_Xenium", pt.size = 0)
#DotPlot(xenium.combined, features = mk)
dev.off()


pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_vln_nFeature_Xenium_p2_cutoff",cut,"_fc_",fc,"_6.pdf"), width = 12, height = 5)
VlnPlot(xenium.combined, features = "nFeature_Xenium", pt.size = 0)
#DotPlot(xenium.combined, features = mk)
dev.off()

#rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.rds"
#data=readRDS(rds)
#marker_genes=c("THSD7B","TENM1","SPRR1B","CYSRT1","ATP8A2","VIM","AQP4","FAM78B","NELL2","EFEMP1","SCGB2A1","CD74","PTPRC","C1QC","TLR2","CEMIP","LUM","IL13RA2","CLMP", "MMP2")

#pdf(paste0(output_PATH, target_tissue, "_oligo_snRNA_dot_p2_cutoff",cut,".pdf"), width = 12, height = 5)
#DotPlot(ref.combined, features = marker_genes, groupby="subclass")
#dev.off()
#ref.combined@assays[["RNA"]]@data
