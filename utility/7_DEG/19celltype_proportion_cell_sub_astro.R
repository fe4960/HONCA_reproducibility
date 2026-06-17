library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
#cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell", "Astrocyte")
cell=c("Astrocyte", "Fibroblast")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")
meta[meta$race=="AsianIndian",]$race="Asian"
df1=NULL
fn_all=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_majorclass_num_atlas_cell_final")
######meta=meta[meta$tissue=="ON",]
data_all=read.table(fn_all,header=F,sep="\t")


for(c in cell){
data_all1=data_all[data_all$V2!=c,]
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)

data=data[data$V1%in%meta$donor,]
cl=unique(data$V2)
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
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3), ct=V2, age=age, race=race, gender=gender)
#df$ct=data$V2
#df$age=data$age
#df$race=data$race
#df$gender=data$gender


#if(c=="Oligodendrocyte"){
df1=df[(df$total>=1000)&(df$ct %in% cl),]
#}

#if(c=="Oligodendrocyte_precursor_cell"){
#df1=df[(df$total>=20)&(df$ct==c),]
	#df=df[df$total>=20,]
#}

df1$majorclass=c
#df1=rbind(df,df1)
mt=matrix(NA,ncol=2,nrow=length(sub))
for(i in sub){
#df1=df[df$ct==i,]
p=summary(lm(prop~age+race+gender,data=df1))
print(p)
}
#geom_bar(stat="identity",position="fill")
########3p=ggplot(data=df1,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df1$ct)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")

p=ggplot(data=df1,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df1$ct)) + theme_minimal(base_size = 30)+scale_color_viridis_d(option = "turbo")+scale_fill_viridis_d(option = "turbo")


##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_lm_oligo100_opc20.pdf")

pdf(out1,width=20)
print(p)
dev.off()
#########
p=ggplot(data=df1,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "loess", se = F, aes(color=ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_loess_oligo100_opc20.pdf")

pdf(out1,width=10)
print(p)
dev.off()

########
df1$age_range=">=30&<70"
df1[df1$age>=70,]$age_range=">=70"
df1[df1$age<30,]$age_range="<30"
df=df1
#df=df1[df1$age_range != "<60",]
#geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2,alpha = 0.7)
#####p=ggplot(data=df1,aes(x=ct,y=prop,color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylim(0,0.3)  #geom_jitter(width=0.2,alpha = 0.7)

#p=ggplot(data=df1,aes(y=prop,x=age_range, color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=c("gold","#74C476", "purple"))+geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.7)+facet_wrap(~ct, ncol=11, scales="free_y")  #geom_jitter(width=0.2,alpha = 0.7)
#df1$ct_short=df1$ct
if(c=="Astrocyte"){
df1$ct_short <- sapply(strsplit(df1$ct, "_"), function(x) {
  paste(x[1], x[2], sep = "_")
})
}

if(c=="Fibroblast"){
df1$ct_short <- sapply(strsplit(df1$ct, "_"), function(x) {
			       paste(x[c(1, length(x))], collapse = "_")
})
}

#p=ggplot(data=df1,aes(y=prop,x=age_range, color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=c("gold","#74C476", "purple"))+geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.7)+facet_wrap(~ct, ncol=11, scales="free_y")  #geom_jitter(width=0.2,alpha = 0.7)

if(c=="Astrocyte" ){
df1$ct_short=factor(df1$ct_short, levels=c("Astro_ON1", "Astro_ON2", "Astro_ON3", "Astro_ON4", "Astro_ON5", "Astro_ONHON1",  "Astro_ONHON2", "Astro_ONH1",  "Astro_ONH2", "Astro_retina1",  "Astro_retina2"))
p=ggplot(data=df1,aes(y=prop,x=age_range, color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=c("gold","#74C476", "purple"))+geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.7)+facet_wrap(~ct_short, ncol=6, scales="free_y")  #geom_jitter(width=0.2,alpha = 0.7)
}else if(c=="Fibroblast"){
p=ggplot(data=df1,aes(y=prop,x=age_range, color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=c("gold","#74C476", "purple"))+geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.7)+facet_wrap(~ct_short, ncol=6, scales="free_y")  #geom_jitter(width=0.2,alpha = 0.7)

}else{

p=ggplot(data=df1,aes(y=prop,x=age_range, color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=c("gold","#74C476", "purple"))+geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.7)+facet_wrap(~ct, ncol=11, scales="free_y")
}

##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_box_oligo100_opc20.pdf")

pdf(out1,width=14, height=8)
print(p)
dev.off()
write.table(df1, file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_box_oligo100_opc20"), sep="\t")
for(ct in unique(df1$ct)){

df10=df1[df1$ct==ct,]
print(ct)
df10$age_range=factor(df10$age_range, levels=c("<30", ">=30&<70", ">=70"))
p=summary(lm(prop~age_range+race+gender, data=df10))


print(p)

df10$age_range=factor(df10$age_range, levels=c(">=30&<70", "<30", ">=70"))
p=summary(lm(prop~age_range+race+gender, data=df10))


print(p)



#p=wilcox.test(df1[(df1$age<30)&(df1$ct==ct),]$prop, df1[(df1$age>=70)&(df1$ct==ct),]$prop,alternative="greater")$p.value
#print(ct)
#print(p)


p=wilcox.test(df[(df$age<30)&(df$ct==ct),]$prop, df[(df$age>=70)&(df$ct==ct),]$prop,alternative="greater")$p.value
print(ct)
print(p)


p=wilcox.test(df[(df$age_range==">=30&<70")&(df$ct==ct),]$prop, df[(df$age>=70)&(df$ct==ct),]$prop,alternative="greater")$p.value
print(ct)
print(p)


p=wilcox.test(df[(df$age_range=="<30")&(df$ct==ct),]$prop, df[(df$age_range==">=30&<70")&(df$ct==ct),]$prop,alternative="greater")$p.value
print(ct)
print(p)

}

}


############
