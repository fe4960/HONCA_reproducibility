library(ComplexHeatmap)
library(viridis)
library(circlize)
args <- commandArgs(trailingOnly = TRUE)
#fn="HCA_ON/data/7_DEG/genexp_donor_subclass_final/Oligodendrocyte_subclass_DEG/Oligodendrocyte/Oligodendrocyte_GO_top15_GO_Biological_Process_2023_age_down_simple.txt"
fn=args[1] #"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_trg500_enrichR/Oligodendrocyte_precursor_cell_GO_top15_MSigDB_Hallmark_2020_up"
data=read.table(fn,sep="\t",header=T,row.names=1)
mat=data.matrix(data) #t(data)
#pdf(paste0("HCA_ON/data/7_DEG/genexp_donor_subclass_final/Oligodendrocyte_subclass_DEG/Oligodendrocyte/Oligodendrocyte_GO_top15_GO_Biological_Process_2023_age_down_simple.pdf"))
#Heatmap(data1)
#dev.off()

# Function to wrap/shorten names (customize width as needed)
shorten_names <- function(names, width = 15) {
  sapply(names, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# Apply the function
short_names <- shorten_names(gsub("\\."," ", colnames(mat)), width = 40)
#pdf(paste0("HCA_ON/data/7_DEG/genexp_donor_subclass_final/Oligodendrocyte_subclass_DEG/Oligodendrocyte/Oligodendrocyte_GO_top15_GO_Biological_Process_2023_age_down_simple.pdf"),height=12,width=10)
pdf(paste0(fn,".pdf"),height=7,width=12)
# Plot with ComplexHeatmap

#col = colorRamp2(
#    c(min(mat),mean(mat), max(mat)),
#    c("darkblue","darkgreen", "yellow")   # You can increase to viridis(5) or viridis(100) for smoother gradients
#  )
mat_scaled <- t(scale(t(mat)))

Heatmap(
#  mat,
   mat_scaled,
#  name = "Z-score",
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  column_labels = short_names,
  col=viridis(100),
  row_names_max_width = unit(10, "cm"),
  column_names_max_height = unit(10, "cm")
)
dev.off()
