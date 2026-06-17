library(Seurat)
#library(ggplot2)
library(dplyr)



args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

#dir.create(output_PATH, recursive = TRUE)

#coembed.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6", ".rds"))

#####coembed.combined=readRDS(paste0(output_PATH, target_tissue, "_8_coembed.combined_inte_harmony_w1", ".rds"))


coembed.combined=readRDS(paste0(output_PATH, target_tissue, "_8_coembed.combined_inte_harmony_w1", ".rds"))

slide_id=unique(coembed.combined@meta.data$slide_id[coembed.combined@meta.data$ident=="Xenium"])

for( id in slide_id){

oligo=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id==id)&(coembed.combined@meta.data$harmony_anno=="Oligodendrocyte")]

#write.table(oligo,file=paste0(output_PATH, target_tissue, "_xenium_",id,"_oligo"))
write.table(oligo,file=paste0(output_PATH, target_tissue, "_xenium_w1_8_",id,"_oligo"))

fibro=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id==id)&(coembed.combined@meta.data$harmony_anno=="Fibroblast")]

#write.table(fibro,file=paste0(output_PATH, target_tissue, "_xenium_",id,"_fibro"))
write.table(fibro,file=paste0(output_PATH, target_tissue, "_xenium_w1_8_",id,"_fibro"))

astro=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id==id)&(coembed.combined@meta.data$harmony_anno=="Astrocyte")]

#write.table(astro,file=paste0(output_PATH, target_tissue, "_xenium_",id,"_astro"))
write.table(astro,file=paste0(output_PATH, target_tissue, "_xenium_w1_8_",id,"_astro"))

astro=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id==id)&(coembed.combined@meta.data$harmony_anno=="Endothelial_cell")]

#write.table(astro,file=paste0(output_PATH, target_tissue, "_xenium_",id,"_endo"))
write.table(astro,file=paste0(output_PATH, target_tissue, "_xenium_w1_8_",id,"_endo"))


}

#oligo1=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_1")&(coembed.combined@meta.data$p5_clusters %in% c(1,12,15,18,23,28,29))]

#oligo2=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_2")&(coembed.combined@meta.data$p5_clusters %in% c(1,12,15,18,23,28,29))]


#opc1=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_1")&(coembed.combined@meta.data$p5_clusters==11)]
#opc2=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_2")&(coembed.combined@meta.data$p5_clusters==11)]

#fibro1=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_1")&(coembed.combined@meta.data$p5_clusters %in% c(4,7,2,17))]
#fibro2=coembed.combined@meta.data$barcode[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id=="PP_2")&(coembed.combined@meta.data$p5_clusters %in% c(4,7,2,17))]


#write.table(oligo1,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_1"))
#write.table(oligo2,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo_2"))

#write.table(opc1,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc_1"))
#write.table(opc2,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc_2"))

#write.table(fibro1,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro_1"))
#write.table(fibro2,file=paste0(output_PATH, target_tissue, "_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro_2"))

