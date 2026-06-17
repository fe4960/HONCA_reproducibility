library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)
#library(arrow)
library(RColorBrewer)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

cut=20 #args[2]

#####process xenium data

label="fibro"  #args[1]

fc=5 #args[3]
#clu="xenium_clusters_p5" #args[4]
clu="xenium_clusters_p8"
sam="BCM_22_0047_ON_ONH_RNA_ds" #args[5]
sc="subclass1"
#sc="subclass1" #args[6]
#label1="seurat_v3_recluster"
label1="seurat_final_6"
#label1="seurat1"
#label1="seurat_v3"
xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_6_xenium_data_list_comb_cutoff20.rds"))


xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))

#xenium.combined <- FindClusters(xenium.combined, resolution = 1.2, cluster.name = "xenium_clusters_12")
#xenium.combined <- FindClusters(xenium.combined, resolution = 1.5, cluster.name = "xenium_clusters_15")
#xenium.combined <- FindClusters(xenium.combined, resolution = 2, cluster.name = "xenium_clusters_2")
#xenium.combined <- FindClusters(xenium.combined, resolution = 0.5, cluster.name = "xenium_clusters_p5")
#saveRDS(xenium.combined, file=paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))

######pdf(paste0(output_PATH, target_tissue,"_",label, "_xenium_umap_cutoff",cut,"_fc_",fc,"_w2.pdf"), width = 6, height = 5)
########DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_12" ) + ggtitle("Xenium Clusters 1.2")
#######DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_15" ) + ggtitle("Xenium Clusters 1.5")
#######DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_2") + ggtitle("Xenium Clusters 2")
#######DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_p8") + ggtitle("Xenium Clusters 0.8")
########DimPlot(xenium.combined, label = T, group.by = "xenium_clusters_1") + ggtitle("Xenium Clusters 1")
########DimPlot(xenium.combined, label = T, group.by = "slide_id") + ggtitle("slide_id")

########dev.off()



clu="xenium_clusters_p8"

#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL
clu1="xenium_clusters_p5"

x4_meta_all=NULL

for(i in 1:5){



xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)



x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

xenium_obj=xenium_obj[,x1@meta.data$barcode]

####RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))
RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc,"_",label1, "_RCTD_results_nonneuron015.rds"))
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc,"_",label1, "_RCTD_results_nonneuron.rds"))

norm_weights = sweep(RCTD@results$weights, 1, rowSums(RCTD@results$weights), '/')
RCTD@results$results_df["max_weight"]=colnames(norm_weights)[max.col(norm_weights, ties.method = "first")]

bc1=rownames(RCTD@results$results_df)
bc2=rownames(xenium_obj@meta.data)
bc3=bc1[bc1%in%bc2]
xenium_obj=xenium_obj[,bc3]
rownames(x1@meta.data)=x1@meta.data$barcode
xenium_obj@meta.data[bc3,clu1]=x1@meta.data[bc3,clu1]

xenium_obj@meta.data[bc3,clu]=x1@meta.data[bc3,clu]

xenium_obj@meta.data[bc3,"rctd1"]=RCTD@results$results_df[bc3,]$first_type
xenium_obj@meta.data[bc3,"rctd2"]=RCTD@results$results_df[bc3,]$second_type
xenium_obj@meta.data[bc3,"spot_class"]=RCTD@results$results_df[bc3,]$spot_class
xenium_obj@meta.data[bc3,"max_weight"]=RCTD@results$results_df[bc3,]$max_weight


#xenium_obj@meta.data[bc3,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
#xenium_obj@meta.data[bc3,"rctd2"]=RCTD@results$results_df[x1@meta.data$barcode,]$second_type
#xenium_obj@meta.data[bc3,"spot_class"]=RCTD@results$results_df[x1@meta.data$barcode,]$spot_class
#xenium_obj@meta.data[bc3,"max_weight"]=RCTD@results$results_df[x1@meta.data$barcode,]$max_weight


#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==4,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==5,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==6,]$xenium_clusters_p2=0


#if( (i==2)){

meta=xenium_obj@meta.data[xenium_obj@meta.data$spot_class=="singlet",]

tb=meta[,c("rctd1",clu1)]
x1_meta_all=rbind(x1_meta_all,tb)


meta1=xenium_obj@meta.data
tb=meta1[,c("max_weight",clu1)]
x4_meta_all=rbind(x4_meta_all,tb)

#}

#tb=table(xenium_obj@meta.data[,c("rctd1",clu)])
#if( (i==2)){
#}

#max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
#max_colnames

#anno=data.frame(clu=rownames(tb),ct=max_colnames)



##########merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_",label, "_rpca_flt",cut,"_fc_",fc,"_",sam,"_",sc,"_slide_",slide_id,".rds"))
#########if(i!=3){
##########xenium_obj@meta.data[,"harmony_clusters_p6"]=merged_obj@meta.data[paste0(slide_id,"_",x1@meta.data$barcode),]$harmony_clusters_p6
###########3xenium_obj@meta.data[,"harmony_anno"]=merged_obj@meta.data[paste0(slide_id,"_",x1@meta.data$barcode),]$harmony_anno
###########}else{
##########sc="subclass1"
#########har_clu="harmony_clusters_2"
##########xenium_obj@meta.data[[har_clu]]=merged_obj@meta.data[paste0(slide_id,"_",x1@meta.data$barcode),]$harmony_clusters_2
###########tb=table(merged_obj@meta.data[,c(har_clu,sc)])
############max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
############anno=data.frame(p2_clu=rownames(tb),ct=max_colnames)
#############xenium_obj@meta.data$harmony_anno="Unk"
#############for( i in anno$p2_clu){
###########xenium_obj@meta.data$harmony_anno[xenium_obj@meta.data$harmony_clusters_2==i]=anno[anno$p2_clu==i,]$ct
#########}
#########}




#tb=table(xenium_obj@meta.data[,c("harmony_clusters_p6",clu)])
#max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
#max_colnames

#anno=data.frame(clu=rownames(tb),ct=max_colnames)

#tb=table(xenium_obj@meta.data[,c("harmony_anno",clu)])

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"



col=c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_",sc,"_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu,"_", label1, "_only_6.pdf"), width = 50, height = 35)

p=ImageDimPlot(object = xenium_obj, group.by = clu, border.size = 0,  axes=TRUE, dark.background=TRUE, cols = "Paired", border.color=NA)
print(p)
dev.off()

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_",sc,"_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu1,"_", label1, "_only_6.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = clu1, border.size = 0,  axes=TRUE, dark.background=TRUE, cols = "Paired", border.color=NA)
print(p)
dev.off()


###### slide 3 ONONH_1
############cropped.coords <- Crop(xenium_obj[["fov"]], y = c(4500, 6500), x = c(2000, 6000), coords = "plot")
#################xenium_obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
##############DefaultBoundary(xenium_obj[["zoom"]]) <- "segmentation"
#ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome"


pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_",sc,"_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu,"_", label1, "_rctd1_only_6.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj,  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "Paired")
#p=ImageDimPlot(object = xenium_obj, group.by = clu, border.size = 0,  axes=TRUE, dark.background=TRUE, cols = "polychrome", border.color=NA)+scale_fill_manual(values=col)

print(p)
dev.off()

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_",sc,"_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu,"_", label1, "_max_weight_6.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj,  group.by = "max_weight", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "Paired")
#p=ImageDimPlot(object = xenium_obj, group.by = clu, border.size = 0,  axes=TRUE, dark.background=TRUE, cols = "polychrome", border.color=NA)+scale_fill_manual(values=col)

print(p)
dev.off()


pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_",sc,"_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu,"_", label1, "_rctd1_only_singlet_6.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj[,xenium_obj@meta.data$spot_class=="singlet"],  group.by = "rctd1", border.size = 0,  axes=TRUE, dark.background=TRUE,  border.color=NA, cols = "Paired")

print(p)
dev.off()


########pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_harmony_anno_only_harmony_clu_15.pdf"), width = 50, height = 35)

#########p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = NA,  axes=TRUE, dark.background=TRUE,  border.color=NA)
##########print(p)
###########dev.off()


}


x2_meta_all=t(table(x1_meta_all))
#write.table(x2_meta_all, file = paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,"_RCTD1", "_", clu, "_", label1, "_2"), sep="\t", quote=F)
write.table(x2_meta_all, file = paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,"_RCTD1", "_", clu1, "_", label1), sep="\t", quote=F)

max_colnames <- colnames(x2_meta_all)[max.col(x2_meta_all, ties.method = "first")]
max_colnames

anno=data.frame(clu=rownames(x2_meta_all),ct=max_colnames)
print(anno)
write.table(anno, file = paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,"_RCTD1", "_", clu1, "_", label1,"_anno"), sep="\t", quote=F)

x2_meta_all=t(table(x4_meta_all))
#write.table(x2_meta_all, file = paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,"_RCTD1_max_weight", "_", clu, "_", label1, "_2"), sep="\t", quote=F)
write.table(x2_meta_all, file = paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,"_RCTD1_max_weight", "_", clu1, "_", label1), sep="\t", quote=F)

max_colnames <- colnames(x2_meta_all)[max.col(x2_meta_all, ties.method = "first")]
max_colnames

anno=data.frame(clu=rownames(x2_meta_all),ct=max_colnames)
print(anno)
write.table(anno, file = paste0(output_PATH, target_tissue,"_",label, "_6_xenium.combined_cutoff",cut,"_fc_",fc,"_RCTD1_max_weight", "_", clu1, "_", label1,"_anno"), sep="\t", quote=F)



#print(max_colnames)


