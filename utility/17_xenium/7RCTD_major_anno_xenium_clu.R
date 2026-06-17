library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)

od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/spacexr/")
dir.create(output_PATH, recursive = TRUE)

xenium_data_list=readRDS(paste0(output_PATH1, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))
label="major"

sc="majorclass"
#HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium.combined_cutoff20.rds
#xenium.combined=readRDS(paste0(output_PATH1, target_tissue, "coembed.combined_inte_harmony", "_fea.rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))
xenium.combined=readRDS(paste0(output_PATH1, "ONONH_comb_xenium.combined_cutoff20.rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))


for(i in 1:length(xenium_data_list)){
#for(i in 1:2){


xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)



x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

#xenium_obj=xenium_obj[,x1@meta.data$barcode]


RCTD=readRDS(paste0(output_PATH1, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))



#id=unique(xenium.combined@meta.data$slide_id)

#x1@meta.data$rctd1=x1@meta.data$majorclass
#bc=x1@meta.data$barcode[x1@meta.data$ident=="Xenium"]
#xenium_obj@meta.data[,clu]=x1@meta.data[,clu]
x1@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
######x2=xenium.combined[,xenium.combined@meta.data$slide=="snRNA"]
#x1@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
#xenium_obj@meta.data[,"rctd2"]=RCTD@results$results_df[x1@meta.data$barcode,]$second_type

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_ident_ct_anno_rctd_slide_",slide_id,".pdf"), width = 5, height = 5)


#pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ident_ct_anno_fea_rctd_slide_",slide_id,".pdf"), width = 5, height = 5)
DimPlot(x1, reduction = "umap", group.by="rctd1", label =TRUE) + NoLegend()
dev.off()

pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_ident_ct_anno_rctd_slide_",slide_id,".pdf"), width = 12, height = 5)

DimPlot(x1, reduction = "umap", group.by="rctd1", label =FALSE) 


#pdf(paste0(output_PATH, target_tissue, "_xenium_umap_facet_ref_xenium_combined_inte_harmony_ident_ct_anno_fea_rctd_slide_snRNA.pdf"), width = 5, height = 5)
#DimPlot(x2, reduction = "umap", group.by="majorclass", split.by="ident", label =TRUE) + NoLegend()
dev.off()

}
