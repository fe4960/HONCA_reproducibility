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

id=na.omit(id)

ni=12
i=1
x1=readRDS( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches",ni,"_clean_major.rds"))

assay="niche"

x_comb=x1[[assay]]@scale.data
seu=x1
for(i in 2:length(id)){
x1=readRDS( paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_", id[i], "_niches",ni,"_clean_major.rds"))
t_comb=x1[[assay]]@scale.data

x_comb=cbind(x_comb,t_comb)

seu=merge(x=seu,y=x1)

}
assay="niche"
data=x_comb
#########p=fviz_nbclust(data, kmeans, method = "wss", k.max = 16, nstart = 25)

########pdf(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_elbow_niches",ni,"_clean_subclass2.pdf"))
########p + ggtitle("Elbow method for K-means")
########dev.off()

# Silhouette
##########p=fviz_nbclust(data, kmeans, method = "silhouette", k.max = 16, nstart = 25)

##########pdf(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_silhouette_niches",ni,"_clean_subclass2.pdf"))
#########p + ggtitle("silhouette method for K-means")
##########dev.off()

############get_best_k <- function(data, k.max = 16, method = "silhouette", nstart = 25) {
########  res <- fviz_nbclust(data, kmeans, method = method, k.max = k.max, nstart = nstart)
#########  df <- res$data
##########  best_k <- df$clusters[which.max(df$y)]
#########  return(best_k)
##########}

##########best_k <- get_best_k(data)
##########cat("Optimal k =", best_k, "\n")
############print(paste0("silhouette\tk=",best_k))

# Gap statistic
########p=fviz_nbclust(x_comb, kmeans, method = "gap_stat", k.max = 16, nstart = 25)

#########pdf(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_gap_stat_niches",ni,"_clean_subclass.pdf"))
########p + ggtitle("gap_stat method for K-means")
########3dev.off()


library(cluster)
library(factoextra)

set.seed(123)

##########gap_stat <- clusGap(
########  data, 
########  FUN = kmeans, 
#######  nstart = 25, 
#######  K.max = 16, 
#######  B = 50   # number of bootstraps
########)

# Visualize


##########best_k <- which.max(gap_stat$Tab[, "gap"])
#########cat("Best number of clusters (Gap Statistic):", best_k, "\n")
#########print(paste0("Best number of clusters (Gap Statistic):", best_k, "\n"))

##########best_k_1se <- maxSE(gap_stat$Tab[,"gap"],
##########                    gap_stat$Tab[,"SE.sim"],
##########                    method = "firstSEmax")
###########cat("Best k (1-SE rule):", best_k_1se, "\n")
############print(paste0("Best k (1-SE rule):", best_k_1se, "\n"))

##########p=fviz_gap_stat(gap_stat)

############pdf(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_gap_stat_niches",ni,"_clean_subclass1.pdf"))
############p + ggtitle("gap_stat method for K-means")
#############dev.off()

#ggsave("elbow_plot.png", p)
#object <- CreateSeuratObject(counts = x_comb)
#assay="niche"
#object[[assay]]@scale.data=x_comb
niches.k=16
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean16"]] <- results[["cluster"]]

niches.k=6
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean6"]] <- results[["cluster"]]

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

niches.k=8
results <- kmeans(
    x = t(data),
    centers = niches.k,
    nstart = 30,
  )

seu[["kmean8"]] <- results[["cluster"]]



#saveRDS(seu, file=paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6_harmony_anno_gap_stat_niches10_clean_merged.rds"))
saveRDS(seu, file=paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1_harmony_anno_gap_stat_niches10_clean_merged_major.rds"))

library(viridis)

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

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean16_clean_major.pdf"), width=5, height=4)
print(p)
dev.off()


dt=table(seu$harmony_anno, seu$kmean6)

colnames(dt)=seq(1,6,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=6"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean6_clean_major.pdf"), width=5, height=3)
print(p)
dev.off()

dt=table(seu$harmony_anno, seu$kmean10)

colnames(dt)=seq(1,10,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=10"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean10_clean_major.pdf"), width=5, height=3)
print(p)
dev.off()

dt=table(seu$harmony_anno, seu$kmean12)

colnames(dt)=seq(1,12,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=12"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean12_clean_major.pdf"), width=5, height=3)
print(p)
dev.off()


dt=table(seu$harmony_anno, seu$kmean8)

colnames(dt)=seq(1,8,1)

p=pheatmap(
  dt,
  color = viridis(100),    # use viridis palette
  scale = "column",          # don't normalize rows or columns
  cluster_rows = T,    # disable row clustering
  cluster_cols = F,    # disable column clustering
  main = paste0("kmean=8"),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

pdf(paste0(output_PATH, target_tissue, "_all_annot_majorclass_heatmap_niche_kmean8_clean_major.pdf"), width=5, height=3)
print(p)
dev.off()

