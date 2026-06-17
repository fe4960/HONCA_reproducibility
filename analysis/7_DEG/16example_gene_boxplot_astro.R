library(ggplot2)
#library(ggforce)
args <- commandArgs(trailingOnly = TRUE)


main_cell="Astrocyte"

#ct=c("OLIGO1_RBFOX1+", "OLIGO2_LRRC7+", "OLIGO1_SVEP1+", "OLIGO1"   )
#ct=c("Fibro_arachnoid_barrier_STXBP6+",
#"Fibro_arachnoid_CLDN11+",
#"Fibro_arachnoid_trabeculae_HMCN1+",
#"Fibro_arachnoid_TRPM3+",
#"Fibro_dura_boarder_SLC47A1+",
#"Fibro_dura_KCNMA1+",
#"Fibro_dura_NOX4+",
#"Fibro_epipial",
#"Fibro_intima_pia",
#"Fibro_RPEchoroid_BMP5+",
#"Fibro_x"
#)
ct=c(
    "Astro_ON1_WNK2+TSHZ2+",
"Astro_ON2_ACTN1+SERPINA3+",
"Astro_ON3_NR4A3+FOS+",
"Astro_ON4_AFF3+DMGDH+",
"Astro_ON5_DPP10+",
"Astro_ONH1_GABBR2+SV2B+",
"Astro_ONH2_PAX5+GABBR2+",
"Astro_ONHON1_SLC4A11+MARCH1+",
"Astro_ONHON2_CST3+APOE+",
"Astro_retina1_PAX5+ME1+",
"Astro_retina2_NLGN1+"
)
od="genexp_donor_subclass_final_new"

df_all=NULL

for( cell in ct){

dirInName=paste0(cell,"_dream")
fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor"
rml=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass/exp_",cell,"_Norm_all"),header=T)
dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/",cell)
indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass"



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

#######metadata$age_range="between40and70"

#######metadata[metadata$age<=40,]$age_range="younger40"

#########metadata[metadata$age>=70,]$age_range="older70"

########metadata$age_range=factor(metadata$age_range, levels=c("younger31","between31and60","older60"))
if(length(metadata[metadata$race=="AsianIndian",]$race) >0){
metadata[metadata$race=="AsianIndian",]$race="Asian"
}
metadata$race=factor(metadata$race,levels=c("White","Black","Asian","Hispanic","AsianIndian"))
metadata$gender=factor(metadata$gender,levels=c("Female","Male"))

rownames(metadata)=metadata$donor

exp=exp[,rownames(metadata)]

#"LAMA1" "ITGA2" "LAMA3" "ITGA6" "LAMA5" "LAMB1" "LAMB2" "LAMB3" "LAMC1”,” LAMC2" "LAMC3"
#ITGA2
#ITGA6
#ITGA9
#ITGAV
#ITGB8
#ITGB4
#CD44
#SV2B
#DAG1
df=data.frame(donor=rep(metadata$donor,11),  age=rep(metadata$age,11), exp=c(as.numeric(exp["PIEZO1",]), as.numeric(exp["PIEZO2",]),  as.numeric(exp["ITGA2",]), as.numeric(exp["ITGA6",]), as.numeric(exp["ITGA9",]), as.numeric(exp["ITGAV",]), as.numeric(exp["ITGB8",]), as.numeric(exp["ITGB4",]), as.numeric(exp["CD44",]), as.numeric(exp["SV2B",]), as.numeric(exp["DAG1",])       ), age=rep(metadata$age,11), gene=c(rep("PIEZO1",length(metadata$age)), rep("PIEZO2",length(metadata$age)), rep("ITGA2",length(metadata$age)), rep("ITGA6",length(metadata$age)),  rep("ITGA9",length(metadata$age)), rep("ITGAV",length(metadata$age)), rep("ITGB8",length(metadata$age)), rep("ITGB4",length(metadata$age)), rep("CD44",length(metadata$age)), rep("SV2B",length(metadata$age)), rep("DAG1",length(metadata$age)) ), cell=rep(cell, length(metadata$age)*11)  )

#df=data.frame(donor=rep(metadata$donor,11), age_range=rep(metadata$age_range,11), age=rep(metadata$age,11), exp=c(as.numeric(exp["PIEZO1",]), as.numeric(exp["PIEZO2",]),  as.numeric(exp["ITGA2",]), as.numeric(exp["ITGA6",]), as.numeric(exp["ITGA9",]), as.numeric(exp["ITGAV",]), as.numeric(exp["ITGB8",]), as.numeric(exp["ITGB4",]), as.numeric(exp["CD44",]), as.numeric(exp["SV2B",]), as.numeric(exp["DAG1",])       ), age=rep(metadata$age,11), gene=c(rep("PIEZO1",length(metadata$age)), rep("PIEZO2",length(metadata$age)), rep("ITGA2",length(metadata$age)), rep("ITGA6",length(metadata$age)),  rep("ITGA9",length(metadata$age)), rep("ITGAV",length(metadata$age)), rep("ITGB8",length(metadata$age)), rep("ITGB4",length(metadata$age)), rep("CD44",length(metadata$age)), rep("SV2B",length(metadata$age)), rep("DAG1",length(metadata$age)) ), cell=rep(cell, length(metadata$age)*11)  )

df_all=rbind(df_all,df)

}



gene="LAMA_receptor"

#p=ggplot(df_all, aes(x = age, y = as.numeric(exp), color = df$gene, group = df$gene)) +
#  geom_point(size = 3, alpha = 0.8) +
#  geom_smooth(method = "loess", se = FALSE, span = 1,   size = 1) +
#  theme_classic(base_size = 25) +
#  scale_color_brewer(palette = "Set2") +
#  labs(color = "Gene", y = "CPM", x = "")+facet_wrap(~cell)
pdf(paste0(indir,"/OLIGO_",gene,"_exp_level_age_loess.pdf"),height=15,width=20)

p=ggplot(df_all, aes(x = age, y = as.numeric(exp), color = gene, group = gene)) +
#  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 1, size = 1) +
  theme_classic(base_size = 25) +
  scale_color_brewer(palette = "Paired") +
  labs(color = "Gene", y = "CPM", x = "") +  facet_wrap(~cell, scales = "free_y")
print(p)
dev.off()

#########pdf(paste0(indir,"/OLIGO_",gene,"_exp_level_age_box.pdf"),height=7,width=12)
#########p=ggplot(df_all, aes(x = gene, y = as.numeric(exp), fill = age)) +
#########  geom_boxplot() +
#########  theme_classic(base_size = 25) +
#########  scale_fill_brewer(palette = "Set2") +
#########  labs(fill = "Age_range", y = "CPM", x = "") +
########  facet_wrap(~cell, scales = "free_y")

##########print(p)
##########dev.off()
##########

