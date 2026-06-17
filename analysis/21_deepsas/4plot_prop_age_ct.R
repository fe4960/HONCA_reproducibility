library(ggplot2)
library(dplyr)
library(ggpubr)

#celltype,donor,age_year,ifsnc
ct="oligo"
#df1=read.table("software/deepsas/outputs/snc_age_ct_prop", header=T, sep=",")
df1=read.table("software/deepsas/outputs/snc_age_ct_prop_sampleid", header=T, sep=",")

p=ggplot(data=df1,aes(x=age_year,y=ifsnc,fill=celltype))+geom_point(aes(color=celltype,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=celltype))+stat_cor(size=8, aes(color=df1$celltype)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("software/deepsas/outputs/age_prop_", ct, "_ggplot_sampleid.pdf")

pdf(out1,width=10)
print(p)
dev.off()


p=ggplot(data=df1,aes(x=celltype,y=ifsnc,fill=celltype))+geom_boxplot(aes(color=celltype))+geom_jitter(width=0.2, size=1, alpha=0.6)+theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("software/deepsas/outputs/age_prop_", ct, "_ggplot_sampleid_box.pdf")

pdf(out1,width=10)
print(p)
dev.off()


df1=read.table("software/deepsas/outputs/snc_age_ct_prop", header=T, sep=",")

p=ggplot(data=df1,aes(x=celltype,y=ifsnc,fill=celltype))+geom_boxplot(aes(color=celltype))+geom_jitter(width=0.2, size=1, alpha=0.6)+theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("software/deepsas/outputs/age_prop_", ct, "_ggplot_donor_box.pdf")

pdf(out1,width=10)
print(p)
dev.off()


df1=read.table("software/deepsas/outputs/snc_age_class_prop_sampleid", header=T, sep=",")

p=ggplot(data=df1,aes(x=age_year,y=ifsnc))+geom_point(aes(size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x )+stat_cor(size=8) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("software/deepsas/outputs/age_prop_", ct, "_class_ggplot_sampleid.pdf")

pdf(out1,width=10)
print(p)
dev.off()


#df1=read.table("software/deepsas/outputs/snc_ct_prop_sampleid", header=T, sep=",")

