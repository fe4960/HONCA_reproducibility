library('edgeR')
#ct=c("Rod", "Cone", "BC", "AC", "HC", "MG", "RGC", "Astrocyte", "Microglia") #,  "RPE")



#ct=c("Oligodendrocyte","Astrocyte","Fibroblast","Rod","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","Melanocyte","HC","Cone","RGC","Schwann_cell")

args=commandArgs(trailingOnly=T)
ct=read.table(args[1],header=F)$V1
seq="all"
#dir="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas/"
dir=args[2] #"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas"

for(cell in ct){
exp=read.table(paste0(dir,"/exp_",cell),header=T)
fListNames=colnames(exp)

fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)

exp1=exp
countMatrix=data.matrix(exp1)
t=quantile(rowMeans(cpm(countMatrix)),seq(0,1,0.05))
######isexpr=rowMeans(cpm(countMatrix)) >= t[15]

######geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = DGEList( countMatrix )

geneExpr = calcNormFactors( geneExpr )

CPM <- cpm(geneExpr, prior.count=0, log=F)


write.table(CPM,file=paste0(dir,"/exp_",cell,"_Norm_",seq),sep="\t",quote=F)


logCPM <- cpm(geneExpr, prior.count=1, log=TRUE)

write.table(logCPM,file=paste0(dir,"/exp_",cell,"_logNorm_",seq),sep="\t",quote=F)


}

