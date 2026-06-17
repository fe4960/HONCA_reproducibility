library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")

df1=NULL
fn_all=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_majorclass_num_atlas_cell_final")
######meta=meta[meta$tissue=="ON",]
data_all=read.table(fn,header=F,sep="\t")


for(c in cell){
data_all1=data_all[data_all$V2!=c,]
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)

data=data[data$V1%in%meta$donor,]
data_all1=data_all1[data_all1$V1%in%meta$donor,]

data=rbind(data,data_all1)
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender

data=data[order(data$age),]
p=ggplot(data=data,aes(x=age,y=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=20) + theme(legend.position="none")
out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_bar.pdf",width=50)
pdf(out)
print(p)
dev.off()
sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3), ct=ct, age=age, race=race, gender=gender)
#df$ct=data$V2
#df$age=data$age
#df$race=data$race
#df$gender=data$gender


#if(c=="Oligodendrocyte"){
df1=df[(df$total>=1000)&(grepl("OLIGO",df$ct)),]
#}

#if(c=="Oligodendrocyte_precursor_cell"){
#df1=df[(df$total>=20)&(df$ct==c),]
	#df=df[df$total>=20,]
#}

df1$majorclass=c
#df1=rbind(df,df1)
mt=matrix(NA,ncol=2,nrow=length(sub))
#for(i in sub){
#df1=df[df$ct==i,]
#summary(lm(prop~age+race+gender,data=df1))
#geom_bar(stat="identity",position="fill")
p=ggplot(data=df1,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df$ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_lm_oligo100_opc20.pdf")

pdf(out1,width=12)
print(p)
dev.off()
#########
p=ggplot(data=df,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "loess", se = F, aes(color=ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_loess_oligo100_opc20.pdf")

pdf(out1,width=12)
print(p)
dev.off()

########
df$age_range="<60"
df[df$age>=70,]$age_range=">=70"
df[df$age<30,]$age_range="<30"
df=df[df$age_range != "<60",]
#geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2,alpha = 0.7)
p=ggplot(data=df,aes(x=ct,y=prop,color=age_range))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2,alpha = 0.7) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_box_oligo100_opc20.pdf")

pdf(out1,width=12)
print(p)
dev.off()

for(ct in unique(df$ct)){

p=wilcox.test(df[(df$age<30)&(df$ct==ct),]$prop, df[(df$age>=70)&(df$ct==ct),]$prop,alternative="greater")$p.value
print(ct)
print(p)

}

}


############
c="Oligodendrocyte"

fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$ct=c

df=data %>% group_by(V1) %>% summarise(total=sum(V3),age=age, race=race,gender=gender,ct=ct)
df=df[df$total>=500,]
#df_sub=df
df=unique(df)
######
c="Oligodendrocyte_precursor_cell"

fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$ct=c

df1=data %>% group_by(V1) %>% summarise(total=sum(V3),age=age, race=race,gender=gender, ct=ct)
df1=df1[df1$total>=100,]
#df_sub=rbind(df_sub,df1)
df1=unique(df1)
########

n=match(df$V1,df1$V1,nomatch=0)

df1_1=df1[n,]

n=match(df1_1$V1,df$V1,nomatch=0)

df_1=df[n,]

df_all=rbind(df1_1,df_1)

df2=df_all %>% group_by(V1) %>% summarise(prop=total/sum(total),total_all=sum(total),age=age, race=race,gender=gender,celltype=ct)





p=ggplot(data=df2,aes(x=age,y=prop,fill=celltype))+geom_point(aes(color=celltype,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=celltype))+stat_cor(size=8, aes(color=df2$celltype)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_oligo_opc_ggplot_lm.pdf")
pdf(out1,width=14)
print(p)
dev.off()
#######

c="Oligodendrocyte"

fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$ct=c

df=data %>% group_by(V1) %>% summarise(total=sum(V3),age=age, race=race,gender=gender,ct=ct,subtype=V2,count=V3)
df=df[df$total>=500,]
#df_sub=df
#df=unique(df)
######
c="Oligodendrocyte_precursor_cell"

fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$ct=c

df1=data %>% group_by(V1) %>% summarise(total=sum(V3),age=age, race=race,gender=gender, ct=ct, subtype=V2,count=V3)
df1=df1[df1$total>=100,]
#df_sub=rbind(df_sub,df1)
#df1=unique(df1)

donor=unique(df$V1)
donor1=unique(df1$V1)

n=match(donor,donor1,nomatch=0)
dn=donor1[n]

df_f=df[df$V1%in%dn,]
df1_f=df1[df1$V1%in%dn,]
df2=rbind(df_f,df1_f)
#########3
df2_c=df2 %>% group_by(V1) %>% summarise(prop=count/sum(count),total_all=sum(total),age=age, race=race,gender=gender,celltype=ct,subtype=subtype) 

p=ggplot(data=df2_c,aes(x=age,y=prop,fill=subtype))+geom_point(aes(color=subtype,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=subtype))+stat_cor(size=8, aes(color=df2_c$subtype)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_oligo_opc_ggplot_lm_subtype.pdf")
pdf(out1,width=12)
print(p)
dev.off()

sub=unique(df2_c$subtype)
mt=matrix(NA,nrow=length(sub),ncol=2)
i=1
for(s in sub){
	df2_c1=df2_c[df2_c$subtype==s,]
	sum=summary(lm(prop~age+race+gender,data=df2_c1))
	p=sum$coefficients[2,4]
	mt[i,1]=s
	mt[i,2]=p
	i=i+1
}
write.table(mt,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_oligo_opc_ggplot_lm_subtype_p",sep="\t",quote=F)

#df2_c1=df2_c[df2_c$celltype=="Oligodendrocyte",]
df2_c1=df2_c[df2_c$subtype%in% c("OLIGO_SVEP1+"),]

p=ggplot(data=df2_c1,aes(x=age,y=prop,fill=subtype))+geom_point(aes(color=subtype,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=subtype))+stat_cor(size=8, aes(color=df2_c1$subtype)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_oligo_ggplot_lm_subtype_sign.pdf")
pdf(out1,width=12)
print(p)
dev.off()

#df2_c1=df2_c[df2_c$celltype=="Oligodendrocyte_precursor_cell",]
df2_c1=df2_c[df2_c$subtype%in% c("OPC","COP_NFOL","COP_LAMA2+"),]

p=ggplot(data=df2_c1,aes(x=age,y=prop,fill=subtype))+geom_point(aes(color=subtype,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=subtype))+stat_cor(size=8, aes(color=df2_c1$subtype)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_opc_ggplot_lm_subtype_sign.pdf")
pdf(out1,width=12)
print(p)
dev.off()


#########
p=ggplot(data=df2,aes(x=age,y=prop,fill=celltype))+geom_bar(stat="identity",position="fill")+ theme_minimal(base_size = 25)+scale_color_brewer(palette="Paired")+scale_fill_brewer(palette="Paired")
#+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=celltype))+stat_cor(size=8, aes(color=df2$celltye)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_oligo_opc_ggplot_bar.pdf")
pdf(out1,width=12)
print(p)
dev.off()
