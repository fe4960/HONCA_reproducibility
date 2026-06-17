 library(ggplot2)
# cell=c("Astrocyte","Endothelial_cell","Fibroblast","Macrophage","Melanocyte","Microglia","Mural_cell","Oligodendrocyte", "Oligodendrocyte_precursor_cell")
# cell=c("RGC")
 ct="RGC"
 df=NULL

 geno=read.table("/dfs3b/ruic20_lab/junw42/human_ret_anc/scripts/HRCAv2_donor_Chen_ds_MMD_uniq.txt", header=T, sep="\t") 

#geno=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid", header=T, sep=",") 
df_list <- list()

#for(ct in cell){
    fn <- paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/data/3_geneExp/retina/major_subclass/exp_",ct,"_Norm_all")

    df1 <- read.table(fn, header=TRUE, sep="\t", comment="")

    df1 <- df1[rownames(df1) == "RNASEH2B", ]

    df1 <- t(df1)   # remove Chr, start, end, TargetID

    df1 <- as.data.frame(df1)
    colnames(df1)="value"
#    df1$donor=rownames(df1)
    df1$donor=gsub("\\.", "-", rownames(df1))
    df1 <- merge(df1, geno, by="donor")

   df1$value <- as.numeric(df1$value) #as.numeric(df1[,1])

   df1$celltype <- ct
#   df1=df1[df1$tissue=="ONH",]
#    df1$slope = data[data$celltype==ct,]$slope
#    df1$pval=data[data$celltype==ct,]$pval_nominal
    df_list[[ct]] <- df1
   d=mean(df1$value)
   print(ct)
   print(d)
   p=lm(value~age+race+gender, data=df1)
   print(summary(p))
#}

df <- do.call(rbind, df_list)

library(dplyr)
#df_all=df
#df_all <- df_all %>%
#    mutate(
#        facet_label = paste0(
#            celltype, " P=", signif(pval,3),

#            "\nslope=", round(slope,3)
#                    )
#    )

#df_all$genotype=gsub("1", "G/A", df_all$genotype)
#df_all$genotype=gsub("0", "G/G", df_all$genotype)
#df_all$genotype=gsub("2", "A/A", df_all$genotype)


library(ggplot2)
library(ggpubr)

p=ggplot(df, aes(x = age, y = value)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) + stat_cor(size=5, label.y=60) + xlab("Age") + ylab("Normalized Expression")+ facet_wrap(~celltype, scale="free_y") + 
    theme_classic() +
    theme(
        legend.position="none",
	axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
    )
#, label.y=0
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/RNASEH2B_DEG_CPM_retina_RGC.pdf", width=4.5, height=4)

#pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/RNASEH2B_DEG.pdf", width=7, height=6)
print(p)
dev.off()

write.table(df, file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/RNASEH2B_DEG_CPM_demo_RGC", sep="\t", quote=F)
