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

label=args[1]

cut=args[2]
fc=args[3]
sc=args[4]
idx=args[5]
#####process xenium data


#clu="xenium_clusters_p3" #args[4]
sam="BCM_22_0047_ON_ONH_RNA_ds" #args[5]




xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))

#xenium.combined=readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

id=unique(xenium.combined@meta.data$slide_id)


id=na.omit(id)


#x_comb=NULL

i=1

slide_id=id[i] #unique(xenium_data_list[[i]]@meta.data$slide_id)

x_comb=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

fname=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron",idx,".rds")

if(idx == "n"){

fname=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds")

#fname=paste0(output_PATH, "/spacexr/",label,"_", slide_id, "_RCTD_results_nonneuron.rds")

	##fname=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds")
}
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron",idx,".rds"))
RCTD=readRDS(fname)
#if(idx == "n"){
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))
#}

x_comb@meta.data[ ,"rctd1"]=RCTD@results$results_df[x_comb@meta.data$barcode,]$first_type
x_comb@meta.data[,"rctd2"]=RCTD@results$results_df[x_comb@meta.data$barcode,]$second_type
x_comb@meta.data[,"spot_class"]=RCTD@results$results_df[x_comb@meta.data$barcode,]$spot_class

x_comb=x_comb[,x_comb@meta.data[,"rctd1"]!="NA"]
x_comb=x_comb[,x_comb@meta.data$spot_class=="singlet"]

for(i in 2:length(id)){

slide_id=id[i] #unique(xenium_data_list[[i]]@meta.data$slide_id)

x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron",idx,".rds"))

fname=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron",idx,".rds")

if(idx == "n"){
#RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))
fname=paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds")
}

RCTD=readRDS(fname)

x1@meta.data[ ,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
x1@meta.data[,"rctd2"]=RCTD@results$results_df[x1@meta.data$barcode,]$second_type
x1@meta.data[,"spot_class"]=RCTD@results$results_df[x1@meta.data$barcode,]$spot_class

x1=x1[,x1@meta.data[,"rctd1"]!="NA"]
x1=x1[,x1@meta.data$spot_class=="singlet"]

x_comb=merge( x=x_comb, y= x1)
}


markers <- FindAllMarkers(x_comb, group.by = "rctd1", only.pos = TRUE)
markers <- markers %>%
group_by(cluster) %>%
dplyr::filter(avg_log2FC > 1)

write.csv(markers, paste0(output_PATH, target_tissue,"_",label,"_",sc,"_",idx,"_markers_rctd.csv"))

markers <- markers  %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 4)

write.csv(markers, paste0(output_PATH, target_tissue,"_",label,"_",sc,"_",idx,"_markers_rctd_top4.csv"))

mk=markers$gene

#avg_exp <- AverageExpression(x_comb, features = mk, factor(mk, levels=mk), return.seurat = FALSE, group.by="rctd1", assays="Xenium", layer = "data")
avg_exp <- AverageExpression(x_comb, features = mk, return.seurat = FALSE, group.by="rctd1", assays="Xenium", layer = "data")

#DoHeatmap( merged_obj2, features=mk, group.by="harmony_anno", assays="Xenium")
#myColor <- viridis::magma(50)
myColor <- viridis::viridis(100)

# Convert to matrix
mat <- avg_exp$Xenium


#mat=mat[mk,]
# Pkatmap
library(pheatmap)

pdf(paste0(output_PATH, target_tissue, "_",label,"_",sc,"_",idx,"_markers_rctd.pdf"))
pheatmap(mat,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         color = myColor,
         fontsize = 10)

pheatmap(mat,
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         color = myColor,
         fontsize = 10)


dev.off()

