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
ct=args[1]
cell=args[2]
od=args[3]
#j=as.numeric(args[3])
dem=args[4]
dirIn <- getwd()
all_data=NULL
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023","GO_Biological_Process_2023","MSigDB_Hallmark_2020")
set.seed(1234567)
#file=paste0("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/",dem,"_DEG_res_cpm_",dem,"_pm_ct_ancFix")

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes
}


##########
#data=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",ct,"_majorclass_DEG/",cell,"/",cell,"_DEG_age"),header=T,sep="\t")
data=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",ct,"_majorclass_sample_DEG/",cell,"/",cell,"_DEG_",dem),header=T,sep="\t")

#data=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",od,"/",ct,"_subclass_DEG/",cell,"/",cell,"_DEG_age"),header=T,sep="\t")
background=rownames(data)
#geneList_up=rownames(data[(data$logFC>0)&(data$adj.P.Val<0.05),])
#geneList_down=rownames(data[(data$logFC<0)&(data$adj.P.Val<0.05),])
###geneList_up=rownames(data[(data$logFC>0)&(data$adj.P.Val<0.05),])
####geneList_down=rownames(data[(data$logFC<0)&(data$adj.P.Val<0.05),])
geneList_up=rownames(data[(data$logFC>0)&(data$qval<0.1),])
geneList_down=rownames(data[(data$logFC<0)&(data$qval<0.1),])


if (websiteLive) {
    enriched2 <- enrichr(geneList_up, dbs, background = background)
}

#dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",ct,"_majorclass_DEG/",cell,"/",cell,"_DEG_GO_",dem,"_up")

dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",ct,"_majorclass_sample_DEG/",cell,"/",cell,"_DEG_GO_",dem,"_up")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


p=plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40,
           y = "Count", orderBy = "Adjusted.P.value")
pdf("MSigDB_Hallmark_2020","_",cell,"_up.pdf")
print(p)
dev.off()


if (websiteLive) {
    enriched2 <- enrichr(geneList_down, dbs, background = background)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",ct,"_majorclass_sample_DEG/",cell,"/",cell,"_DEG_GO_",dem,"_down")
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


p=plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40,
           y = "Count", orderBy = "Adjusted.P.value")
pdf("MSigDB_Hallmark_2020","_",cell,"_down.pdf")
print(p)
dev.off()

