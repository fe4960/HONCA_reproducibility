library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)

# directory containing all enrichment result files
#dir_path <- "HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/*/*DEG_age_cpm1_down/"
#Astro_retina1_PAX5+ME1+/Astro_retina1_PAX5+ME1+_DEG_age_cpm1_down/enrichr_GO_Biological_Process_2023.txt"
files=read.table("HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/DEG_age_cpm1_down_list", header=F)
# get all enrichr GO BP files recursively
#files <- list.files(
#  path = dir_path,
#  pattern = "enrichr_GO_Biological_Process_2023\\.txt$",
#  full.names = TRUE,
#  recursive = TRUE
#)

# function to read one file and extract cell type
read_one_enrichr <- function(f) {
  df <- fread(f, sep = "\t", header = TRUE)
  
  # use parent folder name if files are inside per-celltype folders
  celltype_raw <- basename(dirname(f))
  
  # alternatively, if cell type is in file name, use:
  # celltype_raw <- tools::file_path_sans_ext(basename(f))
  
  # remove DEG suffix and everything after it
  celltype <- str_replace(celltype_raw, "_DEG_age_cpm1_down.*$", "")
  
  df %>% filter(Adjusted.P.value < 0.05) %>%
    select(Term, Combined.Score) %>%
    mutate(celltype = celltype)
}

# combine all files
df_all <- bind_rows(lapply(files$V1, read_one_enrichr))

# if the same Term appears more than once in one cell type, keep the max Combined.Score
df_all2 <- df_all %>%
  group_by(celltype, Term) %>%
  summarise(Combined.Score = max(Combined.Score, na.rm = TRUE), .groups = "drop")

write.table(df_all2, file="HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/GO_DEG_age_cpm1_down", sep="\t", quote=F)

# convert to wide matrix
mat_df <- df_all2 %>%
  pivot_wider(names_from = Term, values_from = Combined.Score)

# convert to matrix
mat <- as.data.frame(mat_df)
rownames(mat) <- mat$celltype
mat$celltype <- NULL
mat <- as.matrix(mat)

# replace missing values with 0
mat[is.na(mat)] <- 0

# optional: keep only terms appearing in at least 2 cell types
# keep_terms <- colSums(mat > 0) >= 2
# mat <- mat[, keep_terms, drop = FALSE]

# optional: keep top variable terms
# term_var <- apply(mat, 2, var)
# mat <- mat[, order(term_var, decreasing = TRUE)[1:min(50, ncol(mat))], drop = FALSE]

# heatmap of raw Combined.Score
pdf("HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/GO_DEG_age_cpm1_down.pdf")
pheatmap(
  mat,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA
)
dev.off()
