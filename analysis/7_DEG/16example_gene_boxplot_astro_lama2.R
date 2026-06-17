library(ggplot2)
#library(ggforce)
args <- commandArgs(trailingOnly = TRUE)

#args <- commandArgs(trailingOnly = TRUE)

cell=args[1]
dirInName=args[2] #"major_flt_cpm05_lib_fix_anc_correct_all_ancFix"  #args[2]
fn=args[3] #"/storage/chentemp/wangj/human_ret_anc/data/HRCA_obs_sample_corrected.txt"
rml=read.table(args[4],header=T)
dirIn=args[5]
indir=args[6]
gene=args[7]
mismatch_sample=c("GSM7553445","GSM7553420","BCM_22_0389_ON_RNA") #,"MMD_22_18633_RNA_macular_NeuN_s2")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
exp=read.table(paste0(indir,"/exp_",cell, "_Norm_all"),header=T)
file_info=read.table(paste0(fn),header=T,sep=",")
file_info=file_info[(file_info$age_year!="Unk")&(file_info$gender!="Unk")&(!(file_info$donor %in% rml$V1))&(!(file_info$donor %in% mismatch_sample)),]

fListNames=colnames(exp)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)

m=match(file_info$donor,fListNames,nomatch=0)

exp=exp[,m]
fListNames=colnames(exp)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)
colnames(exp)=fListNames
m=match(fListNames,file_info$donor,nomatch=0)

metadata=file_info[m,]
metadata$age=as.numeric(metadata$age_year)

metadata$age_range="between31and60"

metadata[metadata$age<31,]$age_range="younger31"

metadata[metadata$age>60,]$age_range="older60"

#metadata=metadata[metadata$age_range!="31-60",]

metadata$age_range=factor(metadata$age_range, levels=c("younger31","between31and60","older60"))

metadata$race=factor(metadata$race,levels=c("White","Black","Asian","Hispanic","AsianIndian"))
metadata$gender=factor(metadata$gender,levels=c("Female","Male"))

#metadata$tissue=factor(metadata$tissue,levels=c("ON","ONH"))

rownames(metadata)=metadata$donor

exp=exp[,rownames(metadata)]


#df=data.frame(donor=metadata$donor, age=metadata$age, piezo2=exp["PIEZO2",])

#####df=data.frame(donor=rep(metadata$donor,5), age=rep(metadata$age,5), piezo2=c(array(exp["PIEZO1",]), array(exp["PIEZO2",]), array(exp["YAZ1",]), array(exp["AKT1",]), array(exp["MAPK1",])  ), age=rep(metadata$age,5), gene=c(rep("PIEZO1",5), rep("PIEZO2",5), rep("YAZ1",5), rep("AKT1",5), rep("MAPK1",5)))

df=data.frame(donor=rep(metadata$donor,4), age=rep(metadata$age,4), piezo2=c(array(exp["PIEZO1",]), array(exp["PIEZO2",]),  array(exp["AKT1",]), array(exp["MAPK1",])  ), age=rep(metadata$age,4), gene=c(rep("PIEZO1",length(metadata$age)), rep("PIEZO2",length(metadata$age)), rep("AKT1",length(metadata$age)), rep("MAPK1",length(metadata$age))))

pdf(paste0(cell,"_",gene,"_exp_level_age_loess.pdf"),height=5,width=5)

#ggplot(df,aes(x=age_range,y=as.numeric(piezo2),color=gene))+geom_point(size=3, alpha=0.8)+geom_line(loss)+theme_classic(base_size = 25)+scale_color_brewer(palette="Set2")+labs(color="Age_range")+ylab("CPM")+xlab("")+theme(legend.position = "none")
p=ggplot(df, aes(x = age, y = as.numeric(piezo2), color = gene, group = gene)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line() +
  geom_smooth(method = "loess", se = FALSE, span = 1, color=gene,  size = 1) +
  theme_classic(base_size = 25) +
  scale_color_brewer(palette = "Set2") +
  labs(color = "Gene", y = "CPM", x = "") +
  theme(legend.position = "none")

print(p)
#ggplot(df,aes(x=age_range,y=as.numeric(piezo2),color=gene))+geom_boxplot()+geom_jitter()+theme_classic(base_size = 25)+scale_color_brewer(palette="Set2")+labs(color="Age_range")+ylab("CPM")+xlab("")+theme(legend.position = "none")
dev.off()
##########

