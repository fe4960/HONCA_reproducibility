library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell", "Astrocyte", "Fibroblast", "Mural_cell", "Endothelial_cell", "Macrophage", "Microglia")
df_all=NULL
for(c in cell){
df1=NULL
#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_majorclass_final/snRNA_donor_num_cell_",c)
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")
meta[meta$race=="AsianIndian",]$race="Asian"
#snRNA_donor_majorclass_num_atlas_cell_final
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_majorclass_num_atlas_cell_final")
######meta=meta[meta$tissue=="ON",]
data=read.table(fn,header=F,sep="\t")
#n=match(meta$sampleid,data$V1,nomatch=0)
data=data[data$V1%in%meta$donor,]
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data=data[order(data$age),]
p=ggplot(data=data,aes(x=age,y=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=20) + theme(legend.position="none")
#out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_bar.pdf",width=50)
out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_",c,"_ggplot_bar.pdf",width=50)
pdf(out)
print(p)
dev.off()
sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3),V3=V3, ct=V2, age=age, race=race, gender=gender)

df1=df[(df$V3>=20)&(df$ct==c),]

df1$age_range <- ">=30&<70"
df1$age_range[df1$age >= 70] <- ">=70"
df1$age_range[df1$age < 30] <- "<30"

x1 <- df1$prop[df1$age < 30]
y1 <- df1$prop[df1$age >= 70]

x1 <- x1[!is.na(x1)]
y1 <- y1[!is.na(y1)]

print(paste0(c, " age < 30 vs. age >= 70"))

p <- wilcox.test(x1, y1, alternative = "greater")$p.value
print(p)

x2 <- df1$prop[df1$age_range == ">=30&<70"]
y2 <- df1$prop[df1$age >= 70]

x2 <- x2[!is.na(x2)]
y2 <- y2[!is.na(y2)]

print(paste0(c, " age >=30&<70 vs. age >=70"))

  p <- wilcox.test(x2, y2, alternative = "greater")$p.value
  print(p)


x2 <- df1$prop[df1$age_range == ">=30&<70"]
y2 <- df1$prop[df1$age < 30]

x2 <- x2[!is.na(x2)]
y2 <- y2[!is.na(y2)]

print(paste0(c, " age < 30 vs. age >=30&<70"))

p <- wilcox.test(x2, y2, alternative = "greater")$p.value
print(p)



mt=matrix(NA,ncol=2,nrow=length(sub))
####df1=df[df$ct==c,]
print(summary(lm(prop~age+gender+race,data=df1)))


p=ggplot(data=df1,aes(x=age,y=prop))+geom_point(aes(size=6))+geom_smooth(method = "lm", se = T, formula = y ~ x, alpha=0.2)+stat_cor(size=8) + ggtitle(c) + theme_minimal(base_size = 30)+xlab("Age")+ylab("Proportion")+theme(legend.position="none")

out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
pdf(out1,width=8)
print(p)
dev.off()


p=ggplot(data=df1,aes(x=age_range,y=prop,fill=age_range))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+scale_fill_manual(values=c("gold","#74C476", "purple"))+labs(x = "Age Range", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), color="black")+ggtitle(c)+theme(legend.position = "none")

write.table(df1, file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_",c,"_ggplot_overall.txt"), sep="\t", quote=F)

out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_",c,"_ggplot_overall.pdf")

pdf(out1,width=6)
print(p)
dev.off()

df_all=rbind( df_all, df1)

}

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/"

p=ggplot(data=df_all,aes(x=age,y=prop))+geom_point(aes(size=6))+geom_smooth(method = "lm", se = T, formula = y ~ x, alpha=0.2)+stat_cor(size=8) +  theme_minimal(base_size = 30)+xlab("Age")+ylab("Cell Proportion")+theme(legend.position="none")+facet_wrap(~ct, ncol=4, scales="free_y" )
        #scale_color_brewer(palette="Set2")#+scale_fill_brewer(palette="Set2")

out1=paste0(dir, "/snRNA_donor_num_cell_all_ggplot_lm.pdf")
pdf(out1,width=20, height=10)
print(p)
dev.off()
#values=c("yellow","#90EE90", "lightbule")
p=ggplot(data=df_all,aes(x=age_range,y=prop,color=age_range))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+labs(x = "Age Range", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), alpha=0.7)+facet_wrap(~ct, ncol=9, scales="free_y")+
  scale_color_manual(
    values = c(
      "<30" = "gold",
      ">=30&<70" = "#74C476",
      ">=70" = "purple"
    )
  )+theme(legend.position = "bottom")

write.table(df1, file=paste0(dir, "/snRNA_donor_num_cell_all_ggplot_overall.txt"), sep="\t", quote=F)

out1=paste0(dir, "/snRNA_donor_num_cell_all_ggplot_overall.pdf")

pdf(out1,width=25, height=7)
print(p)
dev.off()


