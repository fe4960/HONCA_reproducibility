library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)
#library(arrow)

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

########
#####process xenium data

label="oligo"  #args[1]

fc=5 #args[3]
clu="xenium_clusters_1" #args[4]
sam="BCM_22_0047_ON_ONH_RNA_ds" #args[5]
sc="subclass" #args[6]
################
#od="20250620_PP_comb" #args[1]
#target_tissue="ONONH_comb" #args[2]
############
#path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
#wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
#setwd(wd)
#output_PATH <- paste0(wd, "/",target_tissue, "/")
#cut=50


xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))


xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))

#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

for(i in 1:length(id)){



xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)


##############RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id, "_RCTD_results_nonneuron.rds"))
RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id, "_subclass_RCTD_results_nonneuron015.rds"))

x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]
#bc=rownames(RCTD@results$results_df)
xenium_obj=xenium_obj[,x1@meta.data$barcode]
#xenium_obj=xenium_obj[,bc]

xenium_obj@meta.data[,clu]=x1@meta.data[,clu]
xenium_obj@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
xenium_obj@meta.data[,"rctd2"]=RCTD@results$results_df[x1@meta.data$barcode,]$second_type
xenium_obj@meta.data[,"spot_class"]=RCTD@results$results_df[x1@meta.data$barcode,]$spot_class

#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==4,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==5,]$xenium_clusters_p2=3
#xenium_obj@meta.data[xenium_obj@meta.data$xenium_clusters_p2==6,]$xenium_clusters_p2=0

tb=table(xenium_obj@meta.data[,c(clu,"rctd1")])
max_colnames <- colnames(tb)[max.col(tb, ties.method = "first")]
max_colnames

anno=data.frame(clu=rownames(tb),ct=max_colnames)


xenium_obj@meta.data$rctd_anno="Unk"

for( i in anno$clu){
xenium_obj@meta.data$rctd_anno[xenium_obj@meta.data[,clu]==i]=anno[anno$clu==i,]$ct
}

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only_heatmap.pdf"), width = 7, height = 7)
heatmap(tb)
dev.off()

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"





pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_",clu,"_only.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = clu, border.size = 1,  axes=TRUE, dark.background=TRUE, cols = "polychrome", border.color=NA)
print(p)
dev.off()


pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only_majority_vote_anno.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "rctd_anno", border.size = 1,  axes=TRUE, dark.background=TRUE,  border.color=NA) # , cols = "polychrome")
print(p)
dev.off()

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only_subclass.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "rctd1", border.size = 1,  axes=TRUE, dark.background=TRUE, border.color=NA,  cols = "polychrome")
print(p)
dev.off()

xenium_obj=xenium_obj[,xenium_obj@meta.data$spot_class=="singlet"]

pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only_majority_vote_anno_singlet.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "rctd_anno", border.size = 1,  axes=TRUE, dark.background=TRUE,  border.color=NA) # , cols = "polychrome")
print(p)
dev.off()


#pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_harmony_anno_only.pdf"), width = 50, height = 35)
####pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_harmony_anno_only_harmony_clu_15.pdf"), width = 50, height = 35)

#####p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = NA,  axes=TRUE, dark.background=TRUE,  border.color=NA)
#####print(p)
#####dev.off()

xenium_obj@meta.data$y_coor=x1@images$fov$centroids$y


if((i==2) || (i==5)){
xenium_obj@meta.data$y_coor=-xenium_obj@images$fov$centroids$y
}

#library(ggplot2)
############



xenium_obj@meta.data$bin_y <- cut(
  xenium_obj@meta.data$y_coor,
  breaks = 50,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE
)

xenium_obj_meta=xenium_obj@meta.data

################

library(rlang)

library(dplyr)

clu="rctd1"

df <- xenium_obj_meta %>%
  group_by(bin_y, !!sym(clu)  ) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


#ggplot(df,aes(x=bin,y=prop,fill=xenium_clusters_p2))+geom_bar(stat = "identity") +
#  ylab("Proportion") +
#  theme_minimal()

fn=paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",slide_id,"_cutoff",cut,"_fc_",fc,"_rctd1_only_subclass_y_coor_bin.pdf")
p=ggplot(df,aes(x=bin_y,y=prop,fill= !!sym(clu)  ))+geom_bar(stat = "identity") +
  ylab("Proportion") + theme_minimal() + scale_fill_brewer(palette="Set3")
pdf(fn)
print(p)
dev.off()



}

