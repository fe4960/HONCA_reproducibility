library(Seurat)
library(ggplot2)
library(dplyr)



args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")



coembed.combined=readRDS(paste0(output_PATH, target_tissue, "_6_coembed.combined_inte_harmony_w1", ".rds"))

slide_id=unique(coembed.combined@meta.data$slide_id[coembed.combined@meta.data$ident=="Xenium"])

data=NULL

for( id in slide_id){

df=coembed.combined@meta.data[(coembed.combined@meta.data$ident=="Xenium")&(coembed.combined@meta.data$slide_id==id),]

df1=table(df$harmony_anno)

df2=df1[c("Astrocyte","Endothelial_cell", "Fibroblast", "Immune_cell", "Macrophage", "Melanocyte", "Microglia", "Mural_cell", "Oligodendrocyte", "Oligodendrocyte_precursor_cell")]

df3=data.frame(cell_num=df2, cell_class=c("Astrocyte","Endothelial_cell", "Fibroblast", "Immune_cell", "Macrophage", "Melanocyte", "Microglia", "Mural_cell", "Oligodendrocyte", "Oligodendrocyte_precursor_cell"))
df3$slide=id

data=rbind(data,df3)
}

p=ggplot(data=data, aes(y=slide, x=cell_num.Freq, fill=cell_class))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=25) + theme(legend.position="bottom")+scale_fill_brewer(palette="Paired")#+guides(fill = guide_legend(nrow = 1))
	

write.table(data, file=paste0(output_PATH, "cell_prop.txt"), sep="\t", quote=F)

out1=paste0(output_PATH, "cell_prop_6_coembed.pdf")

pdf(out1,height=6, width=12)
print(p)
dev.off()


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

