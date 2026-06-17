library(EWCE)
library(MAGMA.Celltyping)
library(ggplot2)
library(dplyr)
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)
sam=args[1]
label=args[2]
ct=args[3]
#dirIn="/dfs3b/ruic20_lab/junw42/eye_QTL/data/2_GWAS/MAGMA/"
dirIn1="/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/"
dirIn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/13_GWAS/MAGMA"
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
#dir="/dfs3b/ruic20_lab/junw42/human_ret_anc/data/snRNA/"
dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/manuscript/MAGMA/")
rds=paste0(dir,sam,".rds")
obs=read.table(paste0(dir,sam,".obs.gz"),header=T,comment.char="",sep=",",row.names=1)
exp=readRDS(rds)

sce=as.SingleCellExperiment(exp)

#annotLevels <- list(l1=obs$majorclass) #list(l1 = l1, l2 = l2)
annotLevels <- list(l1=obs[,label]) #list(l1 = l1, l2 = l2)

fNames_ALLCELLS <- EWCE::generate_celltype_data(
    exp=sce,
    annotLevels = annotLevels,
    groupName = "retina"
)

retina=load(file=fNames_ALLCELLS)
#file=read.table(paste0(dirIn1,"MAGMA_filelist"))
file=read.table(paste0(dirIn1,"MAGMA_filelist_rm_H36"))

for(i in 1:length(file$V1)){
MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
magma_dirs = file$V1[i],
 ctd = ctd,
  ctd_species = "human", 
  ctd_name = "retina", 
  run_linear = TRUE, 
  run_top10 = TRUE,
  force_new = TRUE)
merged_results <- MAGMA.Celltyping::merge_results(
MAGMA_results = MAGMA_results)
write.table(merged_results,file=paste0(file$V1[i],"MAGMA_",sam),quote=F,sep="\t")
}
