library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(factoextra)
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



id=c("PP_1", "PP_2", "ONH_1", "ONH_2", "ONH_3")
x1_meta_all=NULL

#id=na.omit(id)

ni=10
i=1
x1=readRDS( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_", id[i], "_niches",ni,"_clean_subclass.rds"))

assay="niche"

x_comb=x1[[assay]]@scale.data
seu=x1
for(i in 2:length(id)){
x1=readRDS( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_", id[i], "_niches",ni,"_clean_subclass.rds"))
t_comb=x1[[assay]]@scale.data

x_comb=cbind(x_comb,t_comb)

seu=merge(x=seu,y=x1)

}
assay="niche"
data=x_comb


library(cluster)
library(factoextra)

set.seed(123)

niches.k=12
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean12"]] <- results[["cluster"]]

niches.k=10
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean10"]] <- results[["cluster"]]

niches.k=7
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean7"]] <- results[["cluster"]]


niches.k=16
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean16"]] <- results[["cluster"]]



saveRDS(seu, file=paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_gap_stat_niches10_clean_merged_more.rds"))

library(viridis)


dt=table(seu$harmony_anno, seu$kmean7)

colnames(dt)=seq(1,7,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = T,    # disable column clustering
  main = paste0("kmean=7"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean7_clean_more.pdf"), width=5, height=4)
print(p)
dev.off()


dt=table(seu$harmony_anno, seu$kmean16)

colnames(dt)=seq(1,16,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = T,    # disable column clustering
  main = paste0("kmean=16"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean16_clean_more.pdf"), width=5, height=4)
print(p)
dev.off()


dt=table(seu$harmony_anno, seu$kmean10)

colnames(dt)=seq(1,10,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = T,    # disable column clustering
  main = paste0("kmean=10"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean10_clean_more.pdf"), width=5, height=4)
print(p)
dev.off()

dt=table(seu$harmony_anno, seu$kmean12)

colnames(dt)=seq(1,12,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = T,    # disable column clustering
  main = paste0("kmean=12"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean12_clean_more.pdf"), width=5, height=4)
print(p)
dev.off()

