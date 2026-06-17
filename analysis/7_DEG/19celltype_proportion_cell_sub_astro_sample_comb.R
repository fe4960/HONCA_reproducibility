library(ggplot2)
library(dplyr)
library(ggpubr)
cell=c("Astrocyte")
#meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")

df1=NULL
#fn_all=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_majorclass_num_atlas_cell_final")
#data_all=read.table(fn_all,header=F,sep="\t")


meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid",header=T,sep=",")
fn_all=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_majorclass_num_atlas_cell_final")
data_all=read.table(fn_all,header=F,sep="\t")
meta[meta$race=="AsianIndian",]$race="Asian"

meta$race=factor(meta$race, levels = c("White", "Hispanic", "Asian", "Black"))

tis=c("ONH", "ON")
t="comb"
#for(t in tis){

for(c in cell){
data_all1=data_all[data_all$V2!=c,]
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c)
data=read.table(fn,header=F)

data$V1[data$V1=="BCM_22_0200_ONH_RNA"]="BCM_22_0896_ONH_RNA_s2"
data$V1[data$V1=="BCM_23_0229_ONH_RNA"]="BCM_23_0231_ONH_RNA"
data$V1[data$V1=="BCM_23_0313_ON_RNA"]="MMD_23_20181_ON_RNA"
data$V1[data$V1=="BCM_23_0358_ON_RNA_s2"]="BCM_23_0231_ON_RNA_s2"



data=data[data$V1%in%meta$sampleid,]
cl=unique(data$V2)
data_all1=data_all1[data_all1$V1%in%meta$sampleid,]

data=rbind(data,data_all1)
n=match(data$V1,meta$sampleid,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$tissue=meta[n,]$tissue
#data=data[data$tissue==t,]
data=data[order(data$age),]
p=ggplot(data=data,aes(x=age,y=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=20) + theme(legend.position="none")
out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c,"_comb_ggplot_bar.pdf",width=50)
pdf(out)
print(p)
dev.off()
sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3), ct=V2, age=age, race=race, gender=gender, tissue=tissue)
df1=df[(df$total>=500)&(df$ct %in% cl),]

df1$majorclass=c
#df1=rbind(df,df1)
mt=matrix(NA,ncol=2,nrow=length(sub))
for(i in sub){
#df1=df[df$ct==i,]
df2=df[df$ct==i,]
n1=length(unique(df2$age))
n2=length(unique(df2$race))
n3=length(unique(df2$gender))
if((n1>=2)&(n2>=2)&(n3>=2)){
p=summary(lm(prop~age+race+gender+tissue,data=df2))
print(t)
print(i)
print(p)
}
}
#p=ggplot(data=df1,aes(x=age,y=prop,color=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df1$ct)) + theme_minimal(base_size = 30)+scale_color_viridis_d(option = "turbo")+scale_fill_viridis_d(option = "turbo")+facet_wrap(~ct, ncol=6, scales="free_y")

df1=df1[df1$ct!="Fibro_RPEchoroid_BMP5+",]

p=ggplot(data=df1,aes(x=age,y=prop,color=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df1$ct)) + theme_minimal(base_size = 30)+scale_color_brewer(palette="Paired")+ylab("Cell proportion")+xlab("Age (years)")  #+facet_wrap(~ct, ncol=5, scales="free_y")+theme(legend.position="none")

#out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_lm_oligo100_opc20.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_", t, "_ggplot_lm_oligo100_opc20_1.pdf")

pdf(out1,width=10, height=6)
print(p)
dev.off()
#########
p=ggplot(data=df1,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "loess", se = F, aes(color=ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
#out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_loess_oligo100_opc20.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_",t, "_ggplot_loess_oligo100_opc20_1.pdf")

pdf(out1,width=10)
print(p)
dev.off()

########
df1$age_range=">=30&<70"
df1[df1$age>=70,]$age_range=">=70"
df1[df1$age<30,]$age_range="<30"

df=df1 #df1[df1$age_range != "<60",]



p=ggplot(data=df1,aes(x=age_range,y=prop,color=age_range))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+labs(x = "Age Range", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), alpha=0.7)+facet_wrap(~ct, ncol=5, scales="free_y")+
  scale_color_manual(
    values = c(
      "<30" = "gold",
      ">=30&<70" = "#74C476",
      ">=70" = "purple"
    )
  )+theme(legend.position = "bottom")


#out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_ggplot_box_oligo100_opc20.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_", t, "_ggplot_box_oligo100_opc20_1.pdf")



pdf(out1,width=17, height=10)
print(p)
dev.off()
write.table(df1, file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c,"_", t, "_ggplot_box_oligo100_opc20"), sep="\t")
for(ct in unique(df1$ct)){


df10=df1[df1$ct==ct,]
df10$age_range=factor(df10$age_range, levels=c("<30", ">=30&<70", ">=70"))
p=summary(lm(prop~age_range+race+gender+tissue, data=df10))
print(t)
print(ct)
print(p)


p=wilcox.test(df1[(df1$age<30)&(df1$ct==ct),]$prop, df1[(df1$age>=70)&(df1$ct==ct),]$prop,alternative="greater")$p.value
print(paste0(ct, "age<30 & age >=70"))
print(p)

p=wilcox.test(df1[(df1$age<30)&(df1$ct==ct),]$prop, df1[(df1$age_range==">=30&<70")&(df1$ct==ct),]$prop,alternative="greater")$p.value
print(paste0(ct, "age<30 & age >=30&<70"))
print(p)


p=wilcox.test(df1[(df1$age>=70)&(df1$ct==ct),]$prop, df1[(df1$age_range==">=30&<70")&(df1$ct==ct),]$prop,alternative="greater")$p.value
print(paste0(ct, "age >=30&<70 & age>=70"))
print(p)

}
}
#}

############
