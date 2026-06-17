library(ggplot2)
#library(ggforce)
args <- commandArgs(trailingOnly = TRUE)


main_cell="Fibroblast"

#ct=c("OLIGO1_RBFOX1+", "OLIGO2_LRRC7+", "OLIGO1_SVEP1+", "OLIGO1"   )
ct=c("Fibro_arachnoid_barrier_STXBP6+",
"Fibro_arachnoid_CLDN11+",
"Fibro_arachnoid_trabeculae_HMCN1+",
"Fibro_arachnoid_TRPM3+",
"Fibro_dura_boarder_SLC47A1+",
"Fibro_dura_KCNMA1+",
"Fibro_dura_NOX4+",
"Fibro_epipial",
"Fibro_intima_pia",
"Fibro_RPEchoroid_BMP5+",
"Fibro_x"
)
od="genexp_donor_subclass_final_new"

df_all=NULL

for( cell in ct){

dirInName=paste0(cell,"_dream")
fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor"
rml=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Fibroblast_subclass/exp_",cell,"_Norm_all"),header=T)
dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Fibroblast_subclass_DEG/",cell)
indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Fibroblast_subclass"



mismatch_sample=c("GSM7553445","GSM7553420","BCM_22_0389_ON_RNA") #,"MMD_22_18633_RNA_macular_NeuN_s2")
#dir.create(dirIn,recursive = TRUE)
#setwd(dirIn)
#dirIn <- getwd()
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

metadata$age_range="between12and70"

metadata[metadata$age<=12,]$age_range="younger12"

metadata[metadata$age>=70,]$age_range="older70"

metadata$age_range=factor(metadata$age_range, levels=c("younger31","between31and60","older60"))
metadata[metadata$race=="AsianIndian",]$race="Asian"
metadata$race=factor(metadata$race,levels=c("White","Black","Asian","Hispanic","AsianIndian"))
metadata$gender=factor(metadata$gender,levels=c("Female","Male"))

rownames(metadata)=metadata$donor

exp=exp[,rownames(metadata)]

#"LAMA1" "LAMA2" "LAMA3" "LAMA4" "LAMA5" "LAMB1" "LAMB2" "LAMB3" "LAMC1”,” LAMC2" "LAMC3"

df=data.frame(donor=rep(metadata$donor,5), age_range=rep(metadata$age_range,5), age=rep(metadata$age,5), exp=c(as.numeric(exp["PIEZO1",]), as.numeric(exp["PIEZO2",]),  as.numeric(exp["LAMA2",]), as.numeric(exp["LAMA4",]), as.numeric(exp["LAMC1",])  ),  gene=c(rep("PIEZO1",length(metadata$age)), rep("PIEZO2",length(metadata$age)), rep("LAMA2",length(metadata$age)), rep("LAMA4",length(metadata$age)),  rep("LAMC1",length(metadata$age))), cell=rep(cell, length(metadata$age)*5)  )

df_all=rbind(df_all,df)

}



gene="LAMA"

#p=ggplot(df_all, aes(x = age, y = as.numeric(exp), color = df$gene, group = df$gene)) +
#  geom_point(size = 3, alpha = 0.8) +
#  geom_smooth(method = "loess", se = FALSE, span = 1,   size = 1) +
#  theme_classic(base_size = 25) +
#  scale_color_brewer(palette = "Set2") +
#  labs(color = "Gene", y = "CPM", x = "")+facet_wrap(~cell)
pdf(paste0(indir,"/OLIGO_",gene,"_exp_level_age_loess.pdf"),height=10,width=17)

p=ggplot(df_all, aes(x = age, y = as.numeric(exp), color = gene, group = gene)) +
#  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 1, size = 1) +
  theme_classic(base_size = 25) +
  scale_color_brewer(palette = "Set2") +
  labs(color = "Gene", y = "CPM", x = "") +  facet_wrap(~cell, scales = "free_y")
print(p)
dev.off()

pdf(paste0(indir,"/OLIGO_",gene,"_exp_level_age_box.pdf"),height=7,width=12)
p=ggplot(df_all, aes(x = gene, y = as.numeric(exp), fill = age_range)) +
  geom_boxplot() +
  theme_classic(base_size = 25) +
  scale_fill_brewer(palette = "Set2") +
  labs(fill = "Age_range", y = "CPM", x = "") +
  facet_wrap(~cell, scales = "free_y")

print(p)
dev.off()
##########

