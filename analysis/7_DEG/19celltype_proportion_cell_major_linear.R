library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
#cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell", "Astrocyte", "Fibroblast", "Mural_cell", "Endothelial_cell", "Macrophage", "Microglia")
cell=c("Microglia")
for(c in cell){
df1=NULL
#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_majorclass_final/snRNA_donor_num_cell_",c)
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid",header=T,sep=",")
meta[meta$race=="AsianIndian",]$race="Asian"
meta$race=factor(meta$race, levels = c("White", "Hispanic", "Asian", "Black"))
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_majorclass_num_atlas_cell_final")
######meta=meta[meta$tissue=="ON",]
data=read.table(fn,header=F,sep="\t")
data$V1[data$V1=="BCM_22_0200_ONH_RNA"]="BCM_22_0896_ONH_RNA_s2"
data$V1[data$V1=="BCM_23_0229_ONH_RNA"]="BCM_23_0231_ONH_RNA"
data$V1[data$V1=="BCM_23_0313_ON_RNA"]="MMD_23_20181_ON_RNA"
data$V1[data$V1=="BCM_23_0358_ON_RNA_s2"]="BCM_23_0231_ON_RNA_s2"
#n=match(meta$sampleid,data$V1,nomatch=0)
data=data[data$V1%in%meta$sampleid,]
n=match(data$V1,meta$sampleid,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$tissue=meta[n,]$tissue

data=data[order(data$age),]
p=ggplot(data=data,aes(x=age,y=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=20) + theme(legend.position="none")
#out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_bar.pdf",width=50)
out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_num_cell_",c,"_ggplot_bar.pdf",width=50)
pdf(out)
print(p)
dev.off()
sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3),V3=V3, ct=V2, age=age, race=race, gender=gender, tissue=tissue)
if(c=="Oligodendrocyte"){
df1=df[(df$V3>=200)&(df$ct=="Oligodendrocyte"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Astrocyte"){
df1=df[(df$V3>=200)&(df$ct=="Astrocyte"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Fibroblast"){
df1=df[(df$V3>=200)&(df$ct=="Fibroblast"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}


if(c=="Mural_cell"){
df1=df[(df$V3>=50)&(df$ct=="Mural_cell"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Endothelial_cell"){
df1=df[(df$V3>=50)&(df$ct=="Endothelial_cell"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}


if(c=="Macrophage"){
df1=df[(df$V3>=20)&(df$ct=="Macrophage"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Microglia"){
df1=df[(df$V3>=20)&(df$ct=="Microglia"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}


if(c=="Oligodendrocyte_precursor_cell"){
df1=df[(df$V3>=20)&(df$ct=="Oligodendrocyte_precursor_cell"),]

df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"


p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)

}
mt=matrix(NA,ncol=2,nrow=length(sub))
####df1=df[df$ct==c,]
print(summary(lm(prop~age+gender+race+tissue,data=df1)))

#####print(summary(lm(prop~age+race+gender,data=df1)))
#p=ggplot(data=df1,aes(x=age,y=prop,fill=ct, color=ct))+geom_point(size=6)+geom_smooth(method = "lm", se = T, formula = y ~ x, alpha=0.2)+stat_cor(size=8) + theme_minimal(base_size = 30) + theme(legend.position="none") #+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")

####p=ggplot(data=df1,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = T, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df1$ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")#+scale_fill_brewer(palette="Set2")


p=ggplot(data=df1,aes(x=age,y=prop))+geom_point(aes(size=6))+geom_smooth(method = "lm", se = T, formula = y ~ x, alpha=0.2)+stat_cor(size=8) + ggtitle(c) + theme_minimal(base_size = 30)+xlab("Age")+ylab("Proportion")+theme(legend.position="none")
	#scale_color_brewer(palette="Set2")#+scale_fill_brewer(palette="Set2")

out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_num_cell_",c,"_ggplot_lm.pdf")
pdf(out1,width=8)
print(p)
dev.off()





#p=ggplot(data=df1,aes(x=age_range,y=prop,fill=age_range))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+scale_fill_manual(values=c("darkgreen","yellow"))+labs(x = "Age Range", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), color="black")+ggtitle(c)+theme(legend.position = "none")


p=ggplot(data=df1,aes(x=age_range,y=prop,fill=age_range))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+scale_fill_manual(values=c("yellow","#018571"))+labs(x = "Age Range", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), color="black")+ggtitle(c)+theme(legend.position = "none")

write.table(df1, file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_num_cell_",c,"_ggplot_overall.txt"), sep="\t", quote=F)

out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_num_cell_",c,"_ggplot_overall.pdf")

pdf(out1,width=6)
print(p)
dev.off()


}



