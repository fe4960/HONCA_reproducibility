library(Seurat)
library(pheatmap)
library(dplyr)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
#id="6" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue1, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)


# Define your marker genes
#marker_genes <- c("Dync1i1", "Grm8", "Ctnna2", "St8sia1", "Nrxn1",
#                  "Ptprr", "Pappa", "Ust", "Gm10754", "Hcn1")

#mk=c("BICC1","GFAP","ST18","HS3ST4","F13A1","CERCAM","PTPRB","ABCC9","GNAT1","ARR3","GRM6","ONECUT1","PAX6","GLUL")
#mk=c("ACKR1","SLC4A11","GRM6","ARR3","PTPRB","BICC1","ONECUT1","CD69","F13A1", "MLANA","RLBP1","CD74","RGS5", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ")


#mk=c("GAD1","SLC4A11","GRM6","ARR3","PTPRB","BICC1","ONECUT1","CD69","F13A1", "MLANA","RLBP1","CD74","RGS5", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ")

mk=c("GAD1","SLC4A11","GRM6","ARR3","PTPRB","BICC1","ONECUT1","CD69","F13A1", "MLANA","RGR","CD74","NOTCH3", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ")

# Extract average expression per cluster

merged_obj=readRDS(paste0(output_PATH, target_tissue, "_6_coembed.combined_inte_harmony_w1", ".rds"))

####merged_obj=readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_p6", ".rds"))
merged_obj2=merged_obj[,merged_obj@meta.data$ident=="Xenium"]


markers <- FindAllMarkers(merged_obj2, group.by = "harmony_anno", only.pos = TRUE)
markers <- markers %>%
group_by(cluster) %>%
dplyr::filter(avg_log2FC > 1)
write.csv(markers, paste0(output_PATH, target_tissue,"_6","_markers_harmony_anno_w1.csv"))

###########write.csv(markers, paste0(output_PATH, target_tissue,"_","_markers_harmony_anno_p6.csv"))


avg_exp <- AverageExpression(merged_obj2, features = factor(mk,levels=mk), return.seurat = FALSE, group.by="harmony_anno", assays="Xenium", layer = "data")

#DoHeatmap( merged_obj2, features=mk, group.by="harmony_anno", assays="Xenium")
#myColor <- viridis::magma(50)
myColor <- viridis::viridis(50)

# Convert to matrix
mat <- avg_exp$Xenium

mat1=mat[,c("AC","Astrocyte","BC","Cone","Endothelial-cell","Fibroblast","HC","Immune-cell","Macrophage","Melanocyte","MG","Microglia","Mural-cell","Oligodendrocyte","Oligodendrocyte-precursor-cell","RGC","RPE","Rod", "Schwann-cell")]


####mat1=mat[,c("AC","Astrocyte","BC","Cone","Endothelial-cell","Fibroblast","HC","Immune-cell","Macrophage","Melanocyte","MG","Microglia","Mural-cell","Oligodendrocyte","OPC","RGC","RPE","Rod", "Schwann-cell")]
mat1=mat1[mk,]
# Plot heatmap
pdf(paste0(output_PATH, target_tissue, "_6_coembed.combined_inte_harmony_mk_w1", ".pdf"))

#####pdf(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_mk_p6", ".pdf"))
pheatmap(mat1,
         scale = "row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
	 color = myColor,
         fontsize = 10)
dev.off()
