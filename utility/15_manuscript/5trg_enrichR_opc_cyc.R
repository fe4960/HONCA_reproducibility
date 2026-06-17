#library("DESeq2")
library("ggplot2")
#library( "gplots" )
#library( "RColorBrewer" )
#library(tidyr)
#library('edgeR')
#library(clusterProfiler)
#library(ggrepel)
#library(org.Hs.eg.db)
library(enrichR)
#library(sva)
args <- commandArgs(trailingOnly = TRUE)
#ct=args[1]
#cell=args[2]
#j=as.numeric(args[3])
#dem=args[4]
#dirIn <- getwd()
#all_data=NULL

set.seed(1234567)
#n=500
n=300
websiteLive <- getOption("enrichR.live")

if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes
}


dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023","GO_Biological_Process_2023","MSigDB_Hallmark_2020")


##########
#ct=c("Fibro_arachnoid_PLCB4+","Fibro_arachnoid_STXBP6+","Fibro_dura","Fibro_perivascular","Fibro_pia_ABCA10+","Fibro_pia_HMCN1+_LC?","Fibro_RPEchoroid_BMP5+","Fibro_RPEchoroid_SMOC2+","Fibro_sclera_KCNMA1+","Fibro_sclera_NOX4+","Fibro_x")
     
##############









#ct=c("COP_LAMA2+","COP_NFOL","MFOL","OPC_GLIS3+","OPC")

celltype1="OPC"
    
ct=c("OPC", "OPC_GLIS3+", "MFOL", "COP_NFOL", "NFOL_UTRN+", "OPC_Cycling")
#ct=c("OPC_Cycling")
celltype="Oligodendrocyte_precursor_cell"
for(c in ct){
#file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/", celltype,"/clean/",celltype,"_subclass_seurat.",c,".trg.csv")

file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/", celltype,"/clean/",celltype,"_subclass_seurat_cycling_rmRPE.",c,".trg.csv")

data=read.table(file,header=T,sep=",")
background=data$gene

data=data[(data$logfoldchanges>0)&(data$pvals_adj<0.01),]
data=data[c(1:n),]
geneList=data$gene
######


if (websiteLive) {
    enriched2 <- enrichr(geneList, dbs, background = background,include_overlap = TRUE)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/",celltype,"/clean/",celltype,"_subclass_seurat_trg",n,"_enrichR/",c)
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

###########printEnrich(enriched2)
}
