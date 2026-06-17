library(ggplot2)
library(dplyr)
library(viridis)
library(circlize)
#library(rrvgo)
#library(org.Hs.eg.db)

eng="enrichr_MSigDB_Hallmark_2020"

#dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_fullcell_10kgene_hs_min_gene_200_hvg_10000_scvi_m"
#dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/major_majorclass_DEG/" #Oligodendrocyte/Oligodendrocyte_DEG_GO_age_cpm1/Oligodendrocyte_GO_BP_gse_age_cpm1
#dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/major_majorclass_sample_DEG/" #Oligodendrocyte/Oligodendrocyte_DEG_GO_age_cpm1/Oligodendrocyte_GO_BP_gse_age_cpm1
sam="donor"
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/"
df_up_all=NULL
df_down_all=NULL
go_up=NULL
go_down=NULL
rm_up=read.table("HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap_enrichr_astro_rm_up", header=F, sep="\t")
rm_down=read.table("HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap_enrichr_astro_rm", header=F, sep="\t")

rm_down$V1 <- trimws(rm_down$V1, which = "right")
rm_up$V1 <- trimws(rm_up$V1, which = "right")


lis=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/ct_list_Astrocyte", header=F)
for( i in lis$V1){
  if(i != "Melanocyte"){

    file=paste0(dir1,"/",i,"/",i, "_DEG_age_cpm1_up/", "enrichr_GO_Biological_Process_2023.txt")
if (file.exists(file)) {
#    file=paste0(dir1,"_",i,"_bg/",eng,".txt")
#    df=read.table(file, header=T, sep="\t")

	txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
  txt <- gsub("\\\\r", "", txt)
  txt <- gsub("\\\\n", "\n", txt)
  txt <- gsub("\\\\t", "\t", txt)

  df <- read.delim(text = txt, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    df=df[!(grepl("Alcoho|Skin|Bone|Biomineral|Ossification|Hepatocyte|Auditory|Sensory|Odontogenesis", df$Term)),]


# example enrichR result
df <- read.delim(file)

# extract GO ID from Term column
#df$GOID <- sub(".*\\((GO:[0-9]+)\\)$", "\\1", df$Term)

# score: higher means more important
#scores <- df$Combined.Score
#names(scores) <- df$GOID

# reduce redundancy
#reduced_terms <- reduceSimMatrix(
#  simMatrix = calculateSimMatrix(
#    df$GOID,
#    orgdb = "org.Hs.eg.db",
#    ont = "BP",
#    method = "Rel"
#  ),
#  scores = scores,
#  threshold = 0.7,
#  orgdb = "org.Hs.eg.db"
#)

#reduced_terms

#df=df[df$GOID %in% reduced_terms,]

    df$Term <- sub("\\(.*", "", df$Term)
 df$Term <- trimws(df$Term, which = "right")
    df=df[!(df$Term %in% rm_up$V1),]
    df$Log10P=-log(df$P.value,base=10)
    df$module=paste0(i)
#    df_up=df[df$NES >0,]
#    df_down=df[df$NES <0,]
    df_up=df

#    df_up1=df_up[c(1:10),]
    df_up1 <- df_up[1:min(5, nrow(df_up)), , drop = FALSE]

    #    df_up1 <- df_up[1:min(10, nrow(df_up)), , drop = FALSE]
#    if(length(df_up1[df_up1$Adjusted.P.value<=0.05,1]) >0){
    if (sum(df_up1$Adjusted.P.value <= 0.05, na.rm = TRUE) > 0) {
    df_up1=df_up1[df_up1$Adjusted.P.value<=0.05, drop=F]
    }
    go_up=c(go_up, df_up1$Term)   
    df_up_all=rbind(df_up_all,df_up)
    
#}
}
file=paste0(dir1,"/",i,"/",i, "_DEG_age_cpm1_down/", "enrichr_GO_Biological_Process_2023.txt")
if (file.exists(file)) {
 #   df=read.table(file,header=T,sep="\t")

	txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
  txt <- gsub("\\\\r", "", txt)
  txt <- gsub("\\\\n", "\n", txt)
  txt <- gsub("\\\\t", "\t", txt)

  df <- read.delim(text = txt, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    df=df[!(grepl("Alcoho|Skin|Bone|Biomineral|Ossification|Hepatocyte|Auditory|Sensory|Odontogenesis", df$Term)),]

#####
 # extract GO ID from Term column
#df$GOID <- sub(".*\\((GO:[0-9]+)\\)$", "\\1", df$Term)

# score: higher means more important
#scores <- df$Combined.Score
#names(scores) <- df$GOID

# reduce redundancy
#reduced_terms <- reduceSimMatrix(
#  simMatrix = calculateSimMatrix(
#    df$GOID,
#    orgdb = "org.Hs.eg.db",
#    ont = "BP",
#    method = "Rel"
#  ),
#  scores = scores,
#  threshold = 0.7,
#  orgdb = "org.Hs.eg.db"
#)

#reduced_terms

#df=df[df$GOID %in% reduced_terms,]


    df$Term <- sub("\\(.*", "", df$Term)
     df$Term <- trimws(df$Term, which = "right")
    df=df[!(df$Term %in% rm_down$V1),]

    df$Log10P=-log(df$P.value,base=10)
    df$module=paste0(i)

    df_down=df
    df_down1=df_down[c(1:5),]

 #   df_down1=df_down[c(1:10),]
    df_down1=df_down1[df_down1$Adjusted.P.value<=0.05,]
    go_down=c(go_down, df_down1$Term )   
    df_down_all=rbind(df_down_all,df_down)
}

  }
}

#go=go[go!="Spermatogenesis"]
#    df_up_all[df_up_all$Log10P>=5, ]$Log10P = 5
#    df_down_all[df_down_all$Log10P>=5, ]$Log10P = 5

df_up_all=df_up_all[(df_up_all$Term  %in% go_up),]
df_down_all=df_down_all[(df_down_all$Term  %in% go_down),]

write.table(df_up_all, file=paste0(dir1,"astro_all_enrichr_MSigDB_up_heatmap_top10_",sam, "top_5"), sep="\t", quote=F)
write.table(df_down_all, file=paste0(dir1,"astro_all_enrichr_MSigDB_down_heatmap_top10_",sam, "top_5"), sep="\t", quote=F)

#write.table(df_up_all, file=paste0(dir1,"_all_enrichr_MSigDB_up_heatmap_top5"), sep="\t", quote=F)
#write.table(df_down_all, file=paste0(dir1,"_all_enrichr_MSigDB_down_heatmap_top5"), sep="\t", quote=F)

library(dplyr)
library(tidyr)
library(tibble)

# define your preferred module order
module_order <- lis$V1 #paste0("Module ", sqi)  # or your custom vector

mat <- df_up_all %>%
  select(module, Term , Combined.Score) %>% 
  distinct() %>%   # avoid duplicates
  group_by(module, Term) %>%
  summarise(Combined.Score = max(Combined.Score, na.rm = TRUE), .groups = "drop") %>% 
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module) %>%
  pivot_wider(
    id_cols = module,
    names_from = Term ,
    values_from = Combined.Score,
    values_fill = 0
  ) %>%
  column_to_rownames("module") %>%
  as.matrix()

# row-wise z-score (optional)
mat_scaled <- scale(mat) # t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

#myColor <- viridis::viridis(100)

# Pkatmap
#library(pheatmap)
#mat_scaled=mat_scaled[,go_up]



library(ComplexHeatmap)



shorten_names <- function(names, width = 15) {
  sapply(names, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# Apply the function
short_names <- shorten_names(gsub("\\."," ", colnames(mat_scaled)), width = 35)
out=paste0(dir1,"astro_all_enrichr_MSigDB_up_heatmap_top10_", sam, "_top_5.pdf")


p=Heatmap(
  mat_scaled,
#  name = "Z-score",
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  column_labels = short_names,
  col=viridis(100, option="A"),
  row_names_max_width = unit(10, "cm"),
  column_names_max_height = unit(10, "cm")
)

pdf(out, height=8,width=24)
print(p)
dev.off()


#####
mat <- df_down_all %>%
  select(module, Term , Combined.Score) %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module) %>%
  pivot_wider(
    id_cols = module,
    names_from = Term ,
    values_from = Combined.Score,
    values_fill = 0
  ) %>%
  column_to_rownames("module") %>%
  as.matrix()

# row-wise z-score (optional)
mat_scaled <- scale(mat) #t(scale(t(mat)))
mat_scaled[is.na(mat_scaled)] <- 0

#myColor <- viridis::viridis(100)

# Pkatmap
#library(pheatmap)
#mat_scaled=mat_scaled[,go_down]



library(ComplexHeatmap)



shorten_names <- function(names, width = 15) {
  sapply(names, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# Apply the function
short_names <- shorten_names(gsub("\\."," ", colnames(mat_scaled)), width =35)
out=paste0(dir1,"astro_all_enrichr_MSigDB_down_heatmap_top10_", sam,"_top_5.pdf")


p=Heatmap(
  mat_scaled,
#  name = "Z-score",
  row_names_gp = gpar(fontsize = 18),
  column_names_gp = gpar(fontsize = 18),
  column_labels = short_names,
  col=viridis(100),
  row_names_max_width = unit(10, "cm"),
  column_names_max_height = unit(10, "cm")
)

pdf(out, height=8,width=20)
print(p)
dev.off()


