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
rml=read.table(args[4],header=T)
dirIn=args[5]
indir=args[6]
mismatch_sample=c("GSM7553445","GSM7553420","BCM_22_0389_ON_RNA") #,"MMD_22_18633_RNA_macular_NeuN_s2")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
exp=read.table(paste0(indir,"/exp_",cell),header=T)
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

metadata$age_range="31-60"

metadata[metadata$age<31,]$age_range="younger31"

metadata[metadata$age>60,]$age_range="older60"

metadata=metadata[metadata$age_range!="31-60",]

metadata$age_range=factor(metadata$age_range, levels=c("younger31","older60"))

metadata$race=factor(metadata$race,levels=c("White","Black","Asian","Hispanic","AsianIndian"))
metadata$gender=factor(metadata$gender,levels=c("Female","Male"))

#metadata$tissue=factor(metadata$tissue,levels=c("ON","ONH"))

rownames(metadata)=metadata$donor

exp=exp[,rownames(metadata)]

isexpr=rowMeans(cpm(exp)) >= 5 #t[15]  #rowMeans >=1

countMatrix=exp[isexpr,]
#countMatrix=exp
#geneExpr = DGEList( countMatrix[isexpr,] )

geneExpr = DGEList( countMatrix )
geneExpr = calcNormFactors( geneExpr )


#########

metadata$libsize=floor(log(geneExpr$samples[metadata$donor,]$lib.size,base=10))


#########3

param = SnowParam(4, "SOCK", progressbar=TRUE)


#form <- ~ age  + gender +  race + tissue +  libsize + (1|batch)  ### all_new_all_mac_clean
form <- ~ age_range + libsize + (1|batch) #libsize  # + gender +  race #+ libsize + (1|batch)  ### all_new_all_mac_clean




vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
fitmm = dream( vobjDream, form, metadata)

#fitmm = dream( vobjDream, form, metadata, L)
fitmm = eBayes(fitmm)

head(fitmm$design, 3)

saveRDS(fitmm,file=paste0(cell,"_fitmm_age_range_30_60_cpm5.rds"))

res=topTable( fitmm, coef='age_rangeolder60', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_age_range_30_60_cpm5"),sep="\t",quote=F)

########form <- ~ libsize + (1|age_range)  + (1|gender) +  (1|race)  + (1|batch)

#########vp <- fitExtractVarPartModel(vobjDream, form, metadata)

#########pdf(paste0(cell,"_fitmm_age_range_30_60_cpm1_varPart.pdf"))
##########plotVarPart(sortCols(vp))
##########dev.off()
#res=topTable( fitmm, coef='genderMale', number=Inf )
#res$qval=qvalue(res$P.Value)$qvalues
#write.table(res,file=paste0(cell,"_DEG_gender"),sep="\t",quote=F)


#####res=topTable( fitmm, coef='tissueONH', number=Inf )
#######res$qval=qvalue(res$P.Value)$qvalues
##########write.table(res,file=paste0(cell,"_DEG_tissueONH"),sep="\t",quote=F)


############res=topTable( fitmm, coef='tissuefovea', number=Inf )
############res$qval=qvalue(res$P.Value)$qvalues
############write.table(res,file=paste0(cell,"_DEG_tissueFovea"),sep="\t",quote=F)


#res=topTable( fitmm, coef='raceBlack', number=Inf )
#res$qval=qvalue(res$P.Value)$qvalues
#write.table(res,file=paste0(cell,"_DEG_raceBlack"),sep="\t",quote=F)

########res=topTable( fitmm, coef='raceAsian', number=Inf )
########res$qval=qvalue(res$P.Value)$qvalues
########write.table(res,file=paste0(cell,"_DEG_raceAsian"),sep="\t",quote=F)

#res=topTable( fitmm, coef='raceHispanic', number=Inf )
#res$qval=qvalue(res$P.Value)$qvalues
#write.table(res,file=paste0(cell,"_DEG_raceHispanic"),sep="\t",quote=F)




################form <- ~ age  + gender + race  + tissue  +   libsize +  (age:gender) + (1|batch) ### all_new_all_mac_clean

#################vobjDream = voomWithDreamWeights( geneExpr, form, metadata, BPPARAM=param ) 
##################fitmm = dream( vobjDream, form, metadata )
####################fitmm = eBayes(fitmm)

#head(fitmm$design, 3)

####################saveRDS(fitmm,file=paste0(cell,"_fitmm_inteAgeGender.rds"))


####################res=topTable( fitmm, coef='age:genderMale', number=Inf )
#####################res$qval=qvalue(res$P.Value)$qvalues
#####################write.table(res,file=paste0(cell,"_DEG_ageGender"),sep="\t",quote=F)


