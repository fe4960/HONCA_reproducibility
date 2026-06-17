library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell", "Astrocyte", "Fibroblast", "Mural_cell", "Endothelial_cell", "Macrophage", "Microglia", "Melanocyte")

for(c in cell){
df1=NULL
#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_majorclass_final/snRNA_donor_num_cell_",c)
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")
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
data[data$race=="AsianIndian",]$race="Asian"
data$race=factor(data$race, levels=c("White", "Black", "Hispanic", "Asian"))
sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3),V3=V3, ct=V2, age=age, race=race, gender=gender)
if(c=="Oligodendrocyte"){
df1=df[(df$V3>=50)&(df$ct=="Oligodendrocyte"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Astrocyte"){
df1=df[(df$V3>=100)&(df$ct=="Astrocyte"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Fibroblast"){
df1=df[(df$V3>=50)&(df$ct=="Fibroblast"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}


if(c=="Mural_cell"){
df1=df[(df$V3>=5)&(df$ct=="Mural_cell"),]
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
df1=df[(df$V3>=5)&(df$ct=="Macrophage"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}

if(c=="Microglia"){
df1=df[(df$V3>=5)&(df$ct=="Microglia"),]
df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"
p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)
}


if(c=="Oligodendrocyte_precursor_cell"){
df1=df[(df$V3>=10)&(df$ct=="Oligodendrocyte_precursor_cell"),]

df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"


p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)

}


if(c=="Melanocyte"){
df1=df[(df$V3>0)&(df$ct=="Melanocyte"),]

df1$age_range="<65"
df1[df1$age>=65,]$age_range=">=65"


p=wilcox.test(df1[df1$age<65,]$prop, df1[df1$age>=65,]$prop,alternative="greater")$p.value
print(c)
print(p)

}


mt=matrix(NA,ncol=2,nrow=length(sub))
####df1=df[df$ct==c,]
print(summary(lm(prop~age+gender+race,data=df1)))

}



