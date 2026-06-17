library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
#cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell")
cell=c("Oligodendrocyte")
#meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor",header=T,sep=",")

df1=NULL


meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid",header=T,sep=",")
#fn_all=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_majorclass_num_atlas_cell_final")
#data_all=read.table(fn_all,header=F,sep="\t")
meta[meta$race=="AsianIndian",]$race="Asian"

meta$race=factor(meta$race, levels = c("White", "Hispanic", "Asian", "Black"))

tis=c("ONH", "ON")

#for(t in tis){




for(c in cell){
#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c)
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c)

#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/snRNA_donor_num_cell_",c)
data=read.table(fn,header=F)
#BCM_22_0200_ONH_RNA
#BCM_23_0229_ONH_RNA
#BCM_23_0313_ON_RNA
#BCM_23_0358_ON_RNA_s2
#BCM_22_0458_ON_RNA
#BCM_22_0698_ONH_RNA
#BCM_22_0896_ONH_RNA_s2
#BCM_23_0231_ONH_RNA
#BCM_23_0231_ON_RNA_s2
#BCM_23_1080_ONH_RNA_s1
#BCM_23_1080_ONH_RNA_s2
#GSM7553444
#MMD_23_20181_ON_RNA

data$V1[data$V1=="BCM_22_0200_ONH_RNA"]="BCM_22_0896_ONH_RNA_s2"
data$V1[data$V1=="BCM_23_0229_ONH_RNA"]="BCM_23_0231_ONH_RNA"
data$V1[data$V1=="BCM_23_0313_ON_RNA"]="MMD_23_20181_ON_RNA"
data$V1[data$V1=="BCM_23_0358_ON_RNA_s2"]="BCM_23_0231_ON_RNA_s2"

n=match(data$V1,meta$sampleid,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$tissue=meta[n,]$tissue

data=data[order(data$age),]
#p=ggplot(data=data,aes(x=age,y=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=20) + theme(legend.position="none")
#out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_bar.pdf",width=50)
#pdf(out)
#print(p)
#dev.off()
#sub=unique(data$V2)
df=data %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3), ct=V2, age=age, race=race, gender=gender, tissue=tissue)
#df$ct=data$V2
#df$age=data$age
#df$race=data$race
#df$gender=data$gender
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
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c,"_ggplot_lm_oligo100.pdf")

pdf(out1,width=12)
print(p)
dev.off()
#########
p=ggplot(data=df,aes(x=age,y=prop,fill=ct))+geom_point(aes(color=ct,size=6))+geom_smooth(method = "loess", se = F, aes(color=ct)) + theme_minimal(base_size = 25)+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")
##out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_cell_",c,"_ggplot_lm.pdf")
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c,"_ggplot_loess_oligo100.pdf")

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
out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c,"_ggplot_box_oligo100.pdf")

pdf(out1,width=15)
print(p)
dev.off()

df$age_range=factor(df$age_range, levels=c("<30", ">=30&<70", ">=70"))
res <- df %>%
  group_by(ct) %>%
  group_modify(~{
    fit <- lm(prop ~ age_range + race + gender +  tissue, data = .x)
    coef_tab <- as.data.frame(summary(fit)$coefficients)
    coef_tab$term <- rownames(coef_tab)

    coef_tab %>%
      filter(grepl("^age_range", term)) %>%
      transmute(
        term,
        beta = Estimate,
        se = `Std. Error`,
        t = `t value`,
        pvalue = `Pr(>|t|)`
      )
  }) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(pvalue, method = "BH"))
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c,"_within_lm_age_range_less30.txt")
write.table(res, file = fn, sep="\t", quote=F)

df$age_range=factor(df$age_range, levels=c( ">=30&<70", "<30", ">=70"))
res <- df %>%
  group_by(ct) %>%
  group_modify(~{
    fit <- lm(prop ~ age_range + race + gender +  tissue, data = .x)
    coef_tab <- as.data.frame(summary(fit)$coefficients)
    coef_tab$term <- rownames(coef_tab)

    coef_tab %>%
      filter(grepl("^age_range", term)) %>%
      transmute(
        term,
        beta = Estimate,
        se = `Std. Error`,
        t = `t value`,
        pvalue = `Pr(>|t|)`
      )
  }) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(pvalue, method = "BH"))
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sampleid_subclass_final_new/snRNA_sampleid_num_cell_",c,"_within_lm_age_range_less30-70.txt")
write.table(res, file = fn, sep="\t", quote=F)



for(ct in unique(df$ct)){

df10=df[df$ct==ct,]
df10$age_range=factor(df10$age_range, levels=c("<30", ">=30&<70", ">=70"))
p=summary(lm(prop~age_range+race+gender, data=df10))
print(t)
print(ct)
print(p)

#p=wilcox.test(df[(df$age<30)&(df$ct==ct),]$prop, df[(df$age>=70)&(df$ct==ct),]$prop,alternative="greater")$p.value
#print(ct)
#print(p)


#p=wilcox.test(df[(df$age_range==">=30&<70")&(df$ct==ct),]$prop, df[(df$age>=70)&(df$ct==ct),]$prop,alternative="greater")$p.value
#print(ct)
#print(p)


#p=wilcox.test(df[(df$age_range=="<30")&(df$ct==ct),]$prop, df[(df$age_range==">=30&<70")&(df$ct==ct),]$prop,alternative="greater")$p.value
#print(ct)
#print(p)


}

}


#}
