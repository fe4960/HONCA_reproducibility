# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)
library(harmony)

options(future.globals.maxSize = 50 * 1024^3)
####outdir <- "/dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/analysis/coembed/harmony"
args <- commandArgs(trailingOnly = TRUE)

od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
#id="3" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0)



ref.combined=readRDS(args[1]) #readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table(args[2], header=T,sep=",",row.names=1) #read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)

ref.combined <- AddMetaData(ref.combined, metadata = meta)

###ref.combined = subset(ref.combined, subset = sampleid == args[5])
#ref.combined = subset(ref.combined, subset = sampleid %in% c("BCM_22_0047_ONH_RNA","BCM_22_0047_ON_RNA") )
ref.combined = subset(ref.combined, subset = sampleid %in% c("BCM_22_0047_ONH_RNA","BCM_22_0047_ON_RNA", "BCM_22_0698_ONH_RNA","BCM_22_0784_ON_RNA") )


cut=as.numeric(args[6])
sam=args[5]
#####process xenium data

label=args[4]
fc=as.numeric(args[7])
sc=args[8]
#i=as.numeric(args[9])
t=args[9]


outdir <- paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/imputation/",t, "/")

xenium.combined=readRDS(paste0(outdir, "/", t, "_integrated_merFISH_astro.rds"))

infered.assay = "enhanced"


#xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))
#features <- rownames(xenium.combined)
features = rownames(xenium.combined@assays[[infered.assay]]$data)
#features <- rownames(xenium_data_list[[i]])
features <- intersect(features, rownames(ref.combined))

#print(i)

#for(i in 1:length(xenium_data_list)){
#t=unique(xenium_data_list[[i]]@meta.data$slide_id)

bc=read.table(file=paste0(args[3],"_",t,"_",label),header=T,row.names=1)
#xenium_data_list[[i]]=subset(xenium_data_list[[i]], subset = (barcode %in% bc$x) & (nCount_Xenium >=cut) & (nFeature_Xenium >=fc) )
xenium.combined=subset(xenium.combined, subset = (barcode %in% bc$x) & (nCount_Xenium >=cut) & (nFeature_Xenium >=fc) )


#xenium.combined=xenium_data_list[[i]][features,]
colnames(xenium.combined)=paste0(t,"_",colnames(xenium.combined))
#}





#xenium_data_list <- lapply(X = xenium_data_list, FUN = function(x) {
######xenium.combined <- NormalizeData(xenium.combined)
#####xenium.combined <- ScaleData(xenium.combined, features = features, verbose = FALSE) #, vars.to.regress = "nCount_Xenium")
#####xenium.combined <- RunPCA(xenium.combined, features = features, verbose = FALSE)
#})

######xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "rpca")

#xenium.anchors <- FindIntegrationAnchors(object.list = xenium_data_list, anchor.features = features, reduction = "cca")
##########xenium.combined <- IntegrateData(anchorset = xenium.anchors) #, k.weight = 50)


new_assay <- CreateAssayObject(counts = xenium.combined@assays[[infered.assay]]$data)

xenium.combined@assays[["integrated"]]=new_assay
Key(xenium.combined[["integrated"]])="integrated_"
xenium.combined@assays[["integrated"]]$data <- xenium.combined@assays[["integrated"]]$counts



#DefaultAssay(xenium.combined) <- "integrated"
######
##ref.combined = ref.combined[features,]
###ref.combined <- NormalizeData(ref.combined)
####ref.combined <- ScaleData(ref.combined, features = features, verbose = FALSE, vars.to.regress = "nCount_RNA")


Idents(ref.combined) <- "sampleid"
ref.list <- SplitObject(ref.combined, split.by = "sampleid")

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- NormalizeData(x)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)

    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#features <- intersect(features, rownames(ref.combined))
#features <- c(features, SelectIntegrationFeatures(object.list = ref.list))

features <-  SelectIntegrationFeatures(object.list = ref.list)
features <- c(features, c("SOX2", "SOX9", "VIM"))
features <- unique(features)

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE )
  x <- RunPCA(x, features = features, verbose = FALSE)
})

ref.anchors <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, reduction = "rpca")
ref.combined <- IntegrateData(anchorset = ref.anchors)

DefaultAssay(ref.combined) <- "integrated"
############


#ref.combined <- RunPCA(ref.combined, features = features, verbose = FALSE)
ref.combined@meta.data$slide_id="snRNA"
#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_xenium.combined_cutoff20.rds"))

ref.combined@meta.data$ident <- "Ref"
xenium.combined@meta.data$ident <- "Xenium"
Idents(ref.combined) <- "ident"
Idents(xenium.combined) <- "ident"
#DefaultAssay(ref.combined)="RNA"
DefaultAssay(xenium.combined)="integrated"

######for multiple sample
#ref.combined@assays[["integrated"]]@counts <- ref.combined@assays[["integrated"]]@data
#for one sample

############new_assay <- CreateAssayObject(counts = ref.combined@assays[["RNA"]]@data)

############ref.combined@assays[["integrated"]]=new_assay
###########Key(ref.combined[["integrated"]])="integrated_"
###########ref.combined@assays[["integrated"]]@data <- ref.combined@assays[["integrated"]]@counts

features = rownames(xenium.combined@assays[[infered.assay]]$data)

#features <- rownames(xenium.combined)
features <- intersect(features, rownames(ref.combined))


#xenium.combined <- NormalizeData(xenium.combined)
xenium.combined <- ScaleData(xenium.combined, features = features, verbose = FALSE) #, vars.to.regress = "nCount_Xenium")
xenium.combined <- RunPCA(xenium.combined, features = features, verbose = FALSE)


#########DefaultAssay(ref.combined)="integrated"

#xenium.combined@assays[["integrated"]]$counts <- xenium.combined@assays[["integrated"]]$data

ref.sub <- subset(ref.combined, features = features)
xenium.sub <- subset(xenium.combined, features = features)

#merged_obj <- merge(x = ref.combined[features,], y = xenium.combined[features,])

merged_obj <- merge(x = ref.sub, y = xenium.sub)

# merged_obj <- NormalizeData(merged_obj)


#merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- Seurat::ScaleData(merged_obj, vars.to.regress = "ident" )
#merged_obj <- Seurat::ScaleData(merged_obj, vars.to.regress = "slide_id")

merged_obj <- RunPCA(merged_obj, features = features)



##########merged_obj <- RunHarmony(merged_obj, "ident")
merged_obj <- RunHarmony(merged_obj, group.by.vars="ident", assay.use="integrated") #, max.iter.harmony=100, lambda=0.2)

merged_obj

merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.5, cluster.name = "harmony_clusters_p5")
merged_obj <- FindClusters(merged_obj, resolution = 0.3, cluster.name = "harmony_clusters_p3")
merged_obj <- FindClusters(merged_obj, resolution = 0.6, cluster.name = "harmony_clusters_p6")
merged_obj <- FindClusters(merged_obj, resolution = 1, cluster.name = "harmony_clusters_1")
merged_obj <- FindClusters(merged_obj, resolution = 1.5, cluster.name = "harmony_clusters_15")
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "harmony_clusters_2")

merged_obj <- merged_obj %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

pdf(paste0(output_PATH, target_tissue, "coembed.combined_harmony_umap_",label, "_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 12, height = 5)
####p1 <- DimPlot(merged_obj, reduction = "umap", group.by = "ident")
####p2 <- DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters")
####p1 + p2
###DimPlot(merged_obj, reduction = "umap", group.by = "Atlas_annotation", label = TRUE,repel = TRUE)
####DimPlot(merged_obj, reduction = "umap", group.by = "deconv_annot", label = TRUE, repel = TRUE)
####DimPlot(merged_obj, reduction = "umap", group.by = "final_annot", label = TRUE,repel = TRUE)

######DimPlot(merged_obj, reduction = "umap", group.by = "spot_class", label = TRUE,repel = TRUE)
#DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p3", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p5", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_p6", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_1", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_15", split.by = "ident", label = T)
DimPlot(merged_obj, reduction = "umap", group.by = "harmony_clusters_2", split.by = "ident", label = T)

#DimPlot(merged_obj, reduction = "umap", group.by = "Atlas_annotation", split.by = "ident", label = T)
dev.off()

#table(merged_obj$Atlas_annotation, merged_obj$harmony_clusters)
#table(merged_obj$deconv_annot, merged_obj$harmony_clusters)

saveRDS(merged_obj, paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.rds"))

merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.rds"))

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by=sc,label=T,split.by="ident")+ NoLegend()
dev.off()


merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_1"] <- "PP_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "PP_2"] <- "PP_2"

merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_1"] <- "ONH_1"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_2"] <- "ONH_2"
merged_obj@meta.data$sampleid[merged_obj@meta.data$slide_id == "ONH_3"] <- "ONH_3"



pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_sampleid_",label,"_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 5, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="sampleid")+ NoLegend()
dev.off()

merged_obj@meta.data$slide="snRNA"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_1"]="xenium_1"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "PP_2"]="xenium_2"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_1"] <- "xenium_3"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_2"] <- "xenium_4"
merged_obj@meta.data$slide[merged_obj@meta.data$slide_id == "ONH_3"] <- "xenium_5"


pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_slide_",label,"_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", split.by="slide")+ NoLegend()
dev.off()

table(merged_obj@meta.data[,c(sc,"harmony_clusters_15")])

tb=table(merged_obj@meta.data[,c("harmony_clusters_15",sc)])
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_heatmap_",label,"_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(p15_clu=rownames(tb),ct=max_colnames)

merged_obj@meta.data$harmony_anno="Unk"

for( i in anno$p15_clu){
merged_obj@meta.data$harmony_anno[merged_obj@meta.data$harmony_clusters_15==i]=anno[anno$p15_clu==i,]$ct
}

#####
pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_anno_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by="harmony_anno",label=F,split.by="slide") #+ NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_subclass_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 15, height = 5)
DimPlot(merged_obj, reduction = "umap", group.by=sc,label=F,split.by="slide") #+ NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_count_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 15, height = 5)
#FeaturePlot(merged_obj, reduction = "umap", features="counts",split.by="slide") #+ NoLegend()
FeaturePlot(merged_obj[,merged_obj@meta.data$slide!="snRNA"], reduction = "umap", features="nCount_Xenium",split.by="slide")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_feature_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 15, height = 5)
#FeaturePlot(merged_obj, reduction = "umap", features="counts",split.by="slide") #+ NoLegend()
FeaturePlot(merged_obj[,merged_obj@meta.data$slide!="snRNA"], reduction = "umap", features="nFeature_Xenium",split.by="slide")
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ct_",label,"_cellsize_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.pdf"), width = 15, height = 5)
#FeaturePlot(merged_obj, reduction = "umap", features="counts",split.by="slide") #+ NoLegend()
FeaturePlot(merged_obj[,merged_obj@meta.data$slide!="snRNA"], reduction = "umap", features="cell_size",split.by="slide")
dev.off()

saveRDS(merged_obj, paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",t,"_impute.rds"))

