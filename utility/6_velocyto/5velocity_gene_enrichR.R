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

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes
}


dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023","GO_Biological_Process_2023","MSigDB_Hallmark_2020")



file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_hvg_2000_raw_velocity_gene.csv")

data=read.table(file,header=T,sep=",")

geneList=data$X
######


if (websiteLive) {
    enriched2 <- enrichr(geneList, dbs, background = background,include_overlap = TRUE)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_enrichR")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

##########


file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/oligo_opc_hvg_2000_raw_new_velocity_gene.csv")

data=read.table(file,header=T,sep=",")

geneList=data$X
######


if (websiteLive) {
    enriched2 <- enrichr(geneList, dbs, background = background, include_overlap = TRUE)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/oligo_opc_hvg_2000_raw_new_enrichR")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

##########


##########


file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/veloc/Microglia_subclass_sb_hvg_2000_raw_velocity_gene.csv")

data=read.table(file,header=T,sep=",")

geneList=data$X
######


if (websiteLive) {
    enriched2 <- enrichr(geneList, dbs, background = background, include_overlap = TRUE)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/veloc/Microglia_subclass_sb_hvg_2000_raw_enrichR")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

##########


file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/veloc/Macrophage_subclass_sb_hvg_2000_raw_velocity_gene.csv")

data=read.table(file,header=T,sep=",")

geneList=data$X
######


if (websiteLive) {
    enriched2 <- enrichr(geneList, dbs, background = background, include_overlap = TRUE)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/veloc/Macrophage_subclass_sb_hvg_2000_raw_enrichR")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40, 
           y = "Count", orderBy = "P.value")

