 library(ggplot2)
# data=read.table("eye_QTL/data/5_tensorqtl/retina/ARHGEF_1_16200034_G_A", header=T, sep="\t")
# cell=c("AC", "Astrocyte", "BC", "HC", "MG", "RGC", "Rod" )
 cell=c("Astrocyte","Endothelial_cell","Fibroblast","Macrophage","Melanocyte","Microglia","Mural_cell","Oligodendrocyte", "Oligodendrocyte_precursor_cell")
 df=NULL
# for(ct in cell){
#	fn=paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/data/3_geneExp/retina/major_subclass/exp_",ct,"_inv_norm_sort.bed.gz")
# 	df1=read.table(fn,header=T, comment="", sep="\t")
#	df1=t(df1[df1$TargetID=="ARHGEF19",-c("X.Chr", "start", "end", "TargetID") ])
#	df1[,2]=as.numeric(df1[,2])
#	df=rbind(df,df1)
#  }

# geno <- read.table("sceQTL_exp.raw", header=TRUE)

#geno <- geno[,c("IID","X1.16200034.G.A_G")]
#colnames(geno) <- c("sample","genotype")

geno=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid", header=T, sep=",") 
df_list <- list()

for(ct in cell){
    fn <- paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/major_majorclass_sample/exp_",ct,"_logNorm_all")

    df1 <- read.table(fn, header=TRUE, sep="\t", comment="")

    df1 <- df1[rownames(df1) == "RNASEH2B", ]

    df1 <- t(df1)   # remove Chr, start, end, TargetID

    df1 <- as.data.frame(df1)
    colnames(df1)="value"
    df1$sampleid=rownames(df1)
    df1 <- merge(df1, geno, by="sampleid")

   df1$value <- as.numeric(df1$value) #as.numeric(df1[,1])

    df1$celltype <- ct
   df1=df1[df1$tissue=="ONH",]
#    df1$slope = data[data$celltype==ct,]$slope
#    df1$pval=data[data$celltype==ct,]$pval_nominal
    df_list[[ct]] <- df1
}

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

p=ggplot(df, aes(x = age_year, y = value, color=celltype)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) + stat_cor(size=5, label.y=0) + xlab("Age") + ylab("Normalized Expression")+ facet_wrap(~celltype) + 
    theme_classic() +
    theme(
        legend.position="none",
	axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
    )


pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/RNASEH2B_DEG.pdf", width=7, height=6)
print(p)
dev.off()


