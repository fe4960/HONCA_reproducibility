library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
#cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell")
cell=c("Oligodendrocyte")
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")
meta[meta$race=="AsianIndian",]$race="Asian"
meta$race=factor(meta$race, levels = c("White", "Hispanic", "Asian", "Black"))

df1=NULL
for(c in cell){
#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c)
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)
n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender

data=data[order(data$age),]
#p=ggplot(data=data,aes(x=age,y=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=20) + theme(legend.position="none")
#out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_bar.pdf",width=50)
#pdf(out)
#print(p)
#dev.off()
#sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3))
df$ct=data$V2
df$age=data$age
df$race=data$race
df$gender=data$gender
#df[df$race=="AsianIndian",]$race="Asian"

if(c=="Oligodendrocyte"){
df=df[df$total>=100,]
}

if(c=="Oligodendrocyte_precursor_cell"){
df=df[df$total>=50,]
}
df$majorclass=c
df1=rbind(df,df1)
mt=matrix(NA,ncol=2,nrow=length(sub))
#for(i in sub){
#df1=df[df$ct==i,]
#summary(lm(prop~age+race+gender,data=df1))
#geom_bar(stat="identity",position="fill")
p=ggplot(data=df,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "lm", se = F, formula = y ~ x, aes(color=ct))+stat_cor(size=8, aes(color=df$ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm_oligo100.pdf")

pdf(out1,width=12)
print(p)
dev.off()
#########
p=ggplot(data=df,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "loess", se = F, aes(color=ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_loess_oligo100.pdf")

pdf(out1,width=12)
print(p)
dev.off()

########
df$age_range=">=30&<70"
df[df$age>=70,]$age_range=">=70"
df[df$age<30,]$age_range="<30"
#df=df[df$age_range != "<60",]
#geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2,alpha = 0.7)
p=ggplot(data=df,aes(x=ct,y=prop,color=age_range))+geom_boxplot(outlier.shape = NA) + theme_minimal(base_size = 25)+scale_color_manual(values=c("gold","#74C476", "purple"))+geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.7)+ggtitle(c) #+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_box_oligo100.pdf")

pdf(out1,width=15)
print(p)
dev.off()

for(ct in unique(df$ct)){

df10=df[df$ct==ct,]

#df10$age_range=factor(df10$age_range, levels=c("<30", ">=30&<70", ">=70"))
df10$age_range=factor(df10$age_range, levels=c( ">=30&<70", "<30", ">=70"))
print(ct)

p=summary(lm(prop~age_range+race+gender, data=df10))
print(p)

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



