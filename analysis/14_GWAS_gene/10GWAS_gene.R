library(pheatmap)

df=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/gene_exp_celltype.txt.gz", sep="\t", header=T)

df_filtered <- df[apply(df, 1, max, na.rm = TRUE) < 5, ]


data1=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/gene_exp_majorclass.txt.gz", sep="\t", header=T)

