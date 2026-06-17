library('variancePartition')
library('edgeR')
library('BiocParallel')
#library("DESeq2")
library("ggplot2")
#library( "gplots" )
library( "RColorBrewer" )
#library( "genefilter" )
#library("ggrepel")
library(tidyr)
#library(peer)
library(qvalue)

flt_count = function(ct_exp, cpm_exp,min_ct, min_cpm,min_samples)
{
exp_ct=as.matrix(ct_exp)
exp_ct_raw=as.matrix(ct_exp)
exp_cpm=as.matrix(cpm_exp)
exp_ct[exp_ct < min_ct]=0
exp_ct[exp_ct >= min_ct]=1
exp_cpm[exp_cpm < min_cpm]=0
exp_cpm[exp_cpm >= min_cpm]=1
exp_f = exp_ct_raw[(rowSums(exp_ct)>= (min_samples * ncol(exp_ct)))&(rowSums(exp_cpm)>= (min_samples * ncol(exp_cpm))), ]
return(exp_f)
}



############
args <- commandArgs(trailingOnly = TRUE)

cell=args[1]
dirInName=args[2] #"major_flt_cpm05_lib_fix_anc_correct_all_ancFix"  #args[2]
fn=args[3] #"/storage/chentemp/wangj/human_ret_anc/data/HRCA_obs_sample_corrected.txt"
#rml=read.table(paste0("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas_sample/exp_",cell,"_TMM_cor_flt")) ### final for all except RGC
rml=read.table(args[4],header=T)
dirIn=args[5]
indir=args[6]
##dirIn=paste0("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/",cell,"_",dirInName) #rowMeans >=5
#mix_sample=c("BCM_22_0200_BCM_23_0491_RNA_fovea","Chen_a_10x_Lobe_D017_18_13_NeuN","Chen_a_10x_Lobe_D017_18_13_Nu","Chen_a_10x_Lobe_D018_26_13_NeuN","Chen_a_10x_Lobe_D018_26_13_Nu","GSM7596855","GSM7596860","GSM7596861","GSM7596862","GSM7596863","GSM7596865")
mismatch_sample=c("GSM7553445","GSM7553420","BCM_22_0389_ON_RNA") #,"MMD_22_18633_RNA_macular_NeuN_s2")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
exp=read.table(paste0(indir,"/exp_",cell),header=T)
#cpm=read.table(paste0(indir,"/exp_",cell,"_CPM"),header=T)
#admix=read.table("/storage/chentemp/wangj/human_ret_anc/data/gatk/vcf/hanh_snRNA_1kg_labWGS_anc.prune.fam_5.Q.corrected")
#hanh_list=read.table("/storage/chen/home/jw29/human_ret_anc/scripts/hanh_list",header=T)
#admix$V2=gsub(".filtered","",admix$V2)
#n=match(hanh_list$sampleid,admix$V2,nomatch=0)
#admix$donor=admix$V2
#admix$donor[n]=hanh_list$donorid
file_info=read.table(paste0(fn),header=T,sep=",")
file_info=file_info[(file_info$age_year!="Unk")&(file_info$gender!="Unk")&(!(file_info$sampleid %in% rml$V1))&(!(file_info$sampleid %in% mismatch_sample)),]


#file_info=file_info[(file_info$age!="Unk")&(file_info$gender!="Unk")&(!(file_info$sampleid %in% rml$V1))&(!(file_info$sampleid %in% mix_sample))&(file_info$donor %in% admix$donor)&(!(file_info$sampleid %in% mismatch_sample)),]
#n=match(file_info$donor,admix$donor,nomatch=0)
#file_info$afr=admix$V7[n]
#file_info$adx=admix$V8[n]
#file_info$eur=admix$V9[n]

fListNames=colnames(exp)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)

m=match(file_info$sampleid,fListNames,nomatch=0)

exp=exp[,m]
#cpm=cpm[,m]
fListNames=colnames(exp)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)
colnames(exp)=fListNames
m=match(fListNames,file_info$sampleid,nomatch=0)

metadata=file_info[m,]
metadata$age=as.numeric(metadata$age_year)

metadata$race=factor(metadata$race,levels=c("White","Black","Asian","Hispanic","AsianIndian"))
metadata$gender=factor(metadata$gender,levels=c("Female","Male"))

metadata$tissue=factor(metadata$tissue,levels=c("ON","ONH"))

rownames(metadata)=metadata$sampleid


isexpr=rowMeans(cpm(exp)) > 0.5 #t[15]  #rowMeans >=1

countMatrix=exp[isexpr,]
#countMatrix=exp
#geneExpr = DGEList( countMatrix[isexpr,] )

geneExpr = DGEList( countMatrix )
geneExpr = calcNormFactors( geneExpr )


#########

metadata$libsize=floor(log(geneExpr$samples[metadata$sampleid,]$lib.size,base=10))


#########3

param = SnowParam(4, "SOCK", progressbar=TRUE)


form <- ~ age  + gender +  race + tissue +  libsize + (1|batch)  ### all_new_all_mac_clean


vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
fitmm = dream( vobjDream, form, metadata)

#fitmm = dream( vobjDream, form, metadata, L)
fitmm = eBayes(fitmm)

head(fitmm$design, 3)

saveRDS(fitmm,file=paste0(cell,"_fitmm.rds"))

res=topTable( fitmm, coef='age', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_age"),sep="\t",quote=F)



res=topTable( fitmm, coef='genderMale', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_gender"),sep="\t",quote=F)


res=topTable( fitmm, coef='tissueONH', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_tissueONH"),sep="\t",quote=F)


############res=topTable( fitmm, coef='tissuefovea', number=Inf )
############res$qval=qvalue(res$P.Value)$qvalues
############write.table(res,file=paste0(cell,"_DEG_tissueFovea"),sep="\t",quote=F)


res=topTable( fitmm, coef='raceBlack', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_raceBlack"),sep="\t",quote=F)

########res=topTable( fitmm, coef='raceAsian', number=Inf )
########res$qval=qvalue(res$P.Value)$qvalues
########write.table(res,file=paste0(cell,"_DEG_raceAsian"),sep="\t",quote=F)

res=topTable( fitmm, coef='raceHispanic', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_raceHispanic"),sep="\t",quote=F)




################form <- ~ age  + gender + race  + tissue  +   libsize +  (age:gender) + (1|batch) ### all_new_all_mac_clean

#################vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
##################fitmm = dream( vobjDream, form, metadata )
####################fitmm = eBayes(fitmm)

#head(fitmm$design, 3)

####################saveRDS(fitmm,file=paste0(cell,"_fitmm_inteAgeGender.rds"))


####################res=topTable( fitmm, coef='age:genderMale', number=Inf )
#####################res$qval=qvalue(res$P.Value)$qvalues
#####################write.table(res,file=paste0(cell,"_DEG_ageGender"),sep="\t",quote=F)


