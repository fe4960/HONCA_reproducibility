library(ggplot2)
library(dplyr)
library(ggpubr)
#cell=c("Astrocyte") #,"Oligodendrocyte","Oligodendrocyte_precursor_cell")
cell=c("Oligodendrocyte","Oligodendrocyte_precursor_cell", "Astrocyte", "Fibroblast", "Mural_cell", "Endothelial_cell", "Macrophage", "Microglia", "Melanocyte")

#for(c in cell){
df1=NULL
#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_majorclass_final/snRNA_donor_num_cell_",c)
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid",header=T,sep=",")
#snRNA_donor_majorclass_num_atlas_cell_final
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_sample_majorclass_num_atlas_cell_final")

#fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_majorclass_num_atlas_cell_final")
######meta=meta[meta$tissue=="ON",]
data=read.table(fn,header=F,sep="\t")
#n=match(meta$sampleid,data$V1,nomatch=0)
#data=data[data$V1%in%meta$donor,]
data=data[data$V1%in%meta$sampleid,]
n=match(data$V1,meta$sampleid,nomatch=0)

#n=match(data$V1,meta$donor,nomatch=0)
data$age=meta[n,]$age_year
data$race=meta[n,]$race
data$gender=meta[n,]$gender
data$tissue=meta[n,]$tissue

#data=data[order(data$age),]
data$race1=data$race
data[data$race=="AsianIndian",]$race1="Asian"

data$race1=factor(data$race1, levels = c("white", "Hispanic", "Asian", "Black"))

tis=c("ONH", "ON")
for(t in tis){
data1=data[data$tissue==t,]
p=ggplot(data=data1,aes(y=race1,x=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=25) + theme(legend.position="bottom")+scale_fill_brewer(palette="Paired")+guides(fill = guide_legend(nrow = 3))
out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_race_ggplot_bar_",t,".pdf")
pdf(out, width=10, height=5)
print(p)
dev.off()

p=ggplot(data=data1,aes(y=gender,x=V3,fill=V2))+geom_bar(stat="identity",position="fill")+theme_minimal(base_size=25) + theme(legend.position="bottom")+scale_fill_brewer(palette="Paired")+guides(fill = guide_legend(nrow = 3))
out=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_gender_ggplot_bar_",t,".pdf")
pdf(out, width=10, height=4)
print(p)
dev.off()

df=data1 %>% group_by(V1) %>% summarise(prop=V3/sum(V3),total=sum(V3),V3=V3, ct=V2, age=age, race=race1, gender=gender)

dn=unique(df[(df$ct=="Melanocyte")&(df$prop>=0.3),]$V1)
df=df[(df$ct %in% cell)&(!(df$V1 %in% dn)), ]

if(t == "ON"){
df=df[df$ct!="Melanocyte",]
}

celltype=unique(df$ct)
for(c10 in celltype){
df10=df[df$ct==c10,]
p=summary(lm(prop~gender+race+age, data=df10))
print(paste0(t," ", c10))
print(p)
}

p=ggplot(data=df, aes(x=race, y=prop, color=race))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+scale_color_manual(values=c("seagreen4", "slateblue2", "gold", "lightcoral"))+labs(x = "Ancestry", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), alpha=0.7)+theme(legend.position = "bottom")+facet_wrap(~ct, ncol=10, scales="free_y")

write.table(df, file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_race_ggplot_overall_",t,".txt"), sep="\t", quote=F)

out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_race_ggplot_overall_",t,".pdf")

pdf(out1,width=24)
print(p)
dev.off()

##########
p=ggplot(data=df, aes(x=gender, y=prop, color=gender))+geom_boxplot(outlier.shape = NA)+ theme_minimal(base_size = 25)+scale_color_manual(values=c( "lightcoral", "slateblue2"))+labs(x = "Sex", y = "Cell Proportion")+geom_jitter(position=position_jitterdodge(), alpha=0.7)+theme(legend.position = "bottom")+facet_wrap(~ct, ncol=10, scales="free_y")

write.table(df, file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_gender_ggplot_overall_",t,".txt"), sep="\t", quote=F)

out1=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/snRNA_donor_num_cell_gender_ggplot_overall_",t,".pdf")

pdf(out1,width=15)
print(p)
dev.off()


}
