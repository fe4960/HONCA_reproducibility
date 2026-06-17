library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)

generate_txt=function(outdir, bname){
file=paste0(outdir, "/", bname, "_sample_number.txt")
df_count=read.table(file, header=T, sep="\t")

#####mat <- data %>% pivot_wider(
#####    names_from = sampleid,
#####    values_from = nCount_RNA,
#####    values_fill = 0
#######  ) %>% column_to_rownames("subclass") %>% as.matrix()

######p=pheatmap(mat, 
#####         scale = "row",          # Standardize rows for visual contrast
#####         cluster_rows = TRUE,
#####	 cluster_cols = TRUE,
#####	 fontsize_row = 10,
####  	 fontsize_col = 8,
####  	 angle_col = 45,
#####         main = bname)


library(pheatmap)

#df_count <- read.table(file, header = TRUE, sep = "\t")

# count matrix
mat <- xtabs(nCount_RNA ~ subclass + sampleid, data = df_count)

# convert to proportions within each sample
mat_prop <- prop.table(mat, margin = 1)

# heatmap
p=pheatmap(
  mat_prop,
#  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE
)


write.table(mat_prop, file=paste0(outdir, "/", bname, "_sample_number_pheatmap.txt"), quote=F, sep="\t")

f=paste0(outdir, "/", bname, "_sample_number_pheatmap.pdf")
pdf(f, width=14, height=5)
print(p)
dev.off()

}


outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
bname="Astrocyte_subclass_new5_clean"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
bname="Oligodendrocyte_precursor_cell_subclass_seurat_cycling"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
bname="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location"
generate_txt(outdir, bname)

outdir='/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/'
bname='Endothelial_cell_subclass_sb_seurat_rmRPE'
generate_txt(outdir, bname)


outdir='/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/'
bname="Mural_cell_subclass_rmRPE"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
bname="Oligodendrocyte_subclass_seurat"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
bname="Macrophage_subclass_sb_clean_rmRPE"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
bname="Microglia_subclass_sb_clean"
generate_txt(outdir, bname)


