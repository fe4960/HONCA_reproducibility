library(ggplot2)
library(dplyr)
library(viridis)
library(circlize)
eng="enrichr_MSigDB_Hallmark_2020"

##dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_permut_hs_min_gene_200_hvg_5000_age_m"
#dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_fullcell_10kgene_hs_min_gene_200_hvg_10000_scvi_m"
#dir1="HCA_ON/data/20_hotspot/Astrocyte_hs_min_gene_100_hvg_5000_scvi_m"
dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_fullcell_5kgene_hs_min_gene_250_hvg_5000_scvi_m"

df_all=NULL

go=NULL
sq=c(1,4,5)
for( i in sq){

#for( i in seq(11,12,1)){
    file=paste0(dir1,"_",i,"_bg/",eng,".txt")
    df=read.table(file,header=T,sep="\t")
#    df$Log10P=-log(df$Adjusted.P.value, base=10)

    df$Log10P=-log(df$P.value,base=10)
    df$module=paste0("Module ",i)

    df1=df[c(1:5),]
    df1=df1[df1$Log10P>=2,]
 #   df1=df1[order(df1$Log10P,decreasing=F),]

    go=c(go, df1$Term)   
    df_all=rbind(df_all,df)
}

go=go[go!="Spermatogenesis"]
df_all=df_all[(df_all$Term %in% go),]


library(dplyr)
library(tidyr)
library(tibble)

# define your preferred module order
module_order <- paste0("Module ", sq)  # or your custom vector
mat <- df_all %>%
  select(module, Term, Log10P) %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module) %>%
  pivot_wider(
    id_cols = module,
    names_from = Term,
    values_from = Log10P,
    values_fill = 0
  ) %>%
  column_to_rownames("module") %>%
  as.matrix()

# row-wise z-score (optional)
mat_scaled <- t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

library(ComplexHeatmap)
#myColor <- viridis::viridis(100)

mat=mat[,go]
mat_scaled=mat_scaled[,go]


shorten_names <- function(names, width = 15) {
  sapply(names, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# Apply the function
short_names <- shorten_names(gsub("\\."," ", colnames(mat_scaled)), width = 40)
out=paste0(dir1,"_all_bg_",eng,"_heatmap_update.pdf")


p=Heatmap(
  mat,
#  name = "Z-score",
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  column_labels = short_names,
  col=viridis(100),
  row_names_max_width = unit(10, "cm"),
  column_names_max_height = unit(10, "cm")
)

pdf(out, height=7,width=5)
print(p)
dev.off()

