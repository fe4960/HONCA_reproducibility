library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
#library(arrow)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
############

cut=20

########
#####process xenium data

################
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

#xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6", ".rds"))
xenium.combined=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1", ".rds"))

id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

id=na.omit(id)

#ni=10
ni=12

for(i in 1:length(id)){
x1=xenium.combined[,xenium.combined@meta.data$slide_id==id[i]]

fov="fov"
if(i>1){
fov=paste0("fov.",i)
}

#x1=x1[,!(x1@meta.data$harmony_anno %in% c("Rod","Cone","BC","HC","MG","RGC","AC", "Melanocyte","Schwann_cell"))]
##########csv=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/astro_",id[i],"_subclass1_RCTD_results_nonneuron.txt")
#########oligo=read.table(csv,header=T, sep="\t")
##########idx2=rownames(oligo)
########3adata_index=x1$barcode
########idx2=idx2[idx2 %in% adata_index]

##########n=match(idx2,adata_index,nomatch=0)

###########x1@meta.data[n,"harmony_anno"]=oligo[idx2,"first_type"]

############del_idx=rownames(oligo[oligo$spot_class!="singlet",])

###########csv=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_",id[i],"_subclass_RCTD_results_nonneuron01.txt")


#csv=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_",id[i],"_subclass_RCTD_results_nonneuron015.txt")
##########oligo1=read.table(csv,header=T, sep="\t")
############idx2=rownames(oligo1)
#adata_index=x1$barcode
#########idx2=idx2[idx2 %in% adata_index]
#######n=match(idx2,adata_index,nomatch=0)

##########x1@meta.data[n,"harmony_anno"]=oligo1[idx2,"first_type"]
#########del_idx1=rownames(oligo1[oligo1$spot_class!="singlet",])
#x1=x1[, !(x1$barcode %in% c(del_idx, del_idx1))]

#########x1 <- subset(x1, subset = !(barcode %in% c(del_idx, del_idx1)))

#x1=subset(x1, subset = barcode %in% c(del_idx, del_idx1), invert = TRUE )

#############x1=x1[,!(x1@meta.data$harmony_anno %in% c("Oligodendrocyte","Astrocyte"))]

x1 <- BuildNicheAssay(object = x1,  group.by = "harmony_anno",  niches.k = ni, neighbors.k = 30, fov=fov )


dt=table(x1$harmony_anno, x1$niches)

write.table(dt, file = paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches",ni,"_clean_major.txt"), sep="\t" , quote=F)

saveRDS(x1, file = paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches",ni,"_clean_major.rds"))


DefaultBoundary(x1[[fov]]) <-  "segmentation"

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_niche",ni,"_clean_major.pdf"), width = 50, height = 35)
p= ImageDimPlot(x1, group.by = "niches", border.color=NA, axes=TRUE, border.size = 0, dark.background = T, cols = "polychrome", fov=fov)
print(p)
dev.off()
library(viridis)
n=100

colnames(dt)=seq(1,ni,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("major_", id[i]),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_", id[i], "_heatmap_niche",ni,"_clean_major.pdf"), width=5, height=3)
print(p)
dev.off()


}
