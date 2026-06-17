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
label=args[2]
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


##########Oligodendrocyte_precursor_cell_permut_hs_gene_module_assignment_min_gene_100_age.csv
data=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/",ct,"_hs_gene_module_assignment_min_gene_",label,"_age.csv"),header=T,sep=",")
background=data$Gene
mod=unique(data$Module)
for(m in mod){
if(m != -1){
geneList_up=data[data$Module==m,"Gene"]
#rownames(data[(data$logFC>0)&(data$adj.P.Val<0.05),])
#geneList_down=rownames(data[(data$logFC<0)&(data$adj.P.Val<0.05),])


if (websiteLive) {
    enriched2 <- enrichr(geneList_up, dbs) #, background = background)
}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/",ct,"_hs_min_gene_",label,"_age_m_",m)
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)

printEnrich(enriched2)


p=plotEnrich(enriched2[["MSigDB_Hallmark_2020"]], showTerms = 20, numChar = 40,
           y = "Count", orderBy = "Adjusted.P.value")
pdf(paste0("MSigDB_Hallmark_2020","_",ct,"_",label,"_age_m_",m,".pdf"))
print(p)
dev.off()

}
}
