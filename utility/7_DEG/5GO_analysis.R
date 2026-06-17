library("DESeq2")
library("ggplot2")
library( "gplots" )
library( "RColorBrewer" )
library(tidyr)
library('edgeR')
library(clusterProfiler)
library(ggrepel)
library(org.Hs.eg.db)
#library(sva)
args <- commandArgs(trailingOnly = TRUE)

ct=args[1] #aste0("Schwann_M_treatment",args[1])
cell=args[2]
#cell=args[2]
#seq=args[2]
#dem=args[3]
#j=as.numeric(args[3])
#dem=args[4]
#dirIn=paste0("/dfs3b/ruic20_lab/junw42/DRG/data/3_DEG/DEG_",ct)
########dirIn=paste0("/dfs3b/ruic20_lab/junw42/DRG/data/3_DEG/genexp_sample_majorclass_final/major_subclass_DEG/", ct)
#dirIn=paste0("/dfs3b/ruic20_lab/junw42/DRG/data/3_DEG/genexp_sample_bulk_final/major_subclass_DEG/", ct)
#Oligodendrocyte_subclass_DEG
dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/",ct,"_subclass_DEG/",cell,"_DEG_GO_age_cpm1")

#dirIn=paste0("/dfs3b/ruic20_lab/junw42/DRG/data/3_DEG/genexp_sample_cellclass/major_subclass_DEG/", ct)
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
all_data=NULL

set.seed(1234567)
#file=paste0("chip_",ct,"_DEGs_bulk.txt")
#Macrophage_DEG_surgery_group3

#t=c("chip", "group3")
#for(tt in t){
#file=paste0(cell,"_DEG_surgery_", tt)
fl=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/",ct,"_subclass_DEG/", cell, "/", cell, "_DEG_age_cpm1")
data=read.table(fl, header=T, sep="\t")

#datalist=unique(data$gRNA) #rownames(data)

datalist=rownames(data)

engene_entrez= bitr(datalist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#new_engene=data.frame(SYMBOL=engene_entrez$SYMBOL, ENTREZID=engene_entrez$ENTREZID)

data=data[rownames(data) %in% engene_entrez$SYMBOL,]

#n=match(engene_entrez$SYMBOL,datalist,nomatch=0)


#geneList=data[n,]$log2FoldChange

#geneList=data[n,]$logFC

geneList=data$logFC


m=match(rownames(data),engene_entrez$SYMBOL,nomatch=0)


names(geneList)=engene_entrez[m,]$ENTREZID


#geneList$name=names(geneList)

#geneList1=data.frame(name=engene_entrez[m,]$ENTREZID, fc=geneList)

#geneList1=unique(geneList)


unique_geneList <- geneList[!duplicated(names(geneList))]
geneList=sort(unique_geneList, decreasing = TRUE)

#geneList=geneList[names(geneList)!="NA"]

geneList=geneList[!is.na(names(geneList))]
#############
#GO analysis
###########
library(data.table)
library('biomaRt')
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)


library(GOSemSim)
library(msigdbr)
set.seed(12345)

c="C5"
m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = c , subcategory="GO:BP") %>%
  dplyr::select(gs_name, entrez_gene)

em2 <- GSEA(geneList, TERM2GENE = m_t2g,minGSSize=10,  seed=TRUE, eps=0,nPermSimple = 100000)

if(dim(as.data.frame(em2))[1]!=0){

em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","msigdbr_gse"),quote=F,sep="\t")
#em2_new=read.table(paste0(cell,"_","msigdbr_gse"),header=T,sep="\t")
pdf(paste0(cell,"_","gsea_msigdbr_dotplot","_top10.pdf"),width=8, height=9)
p=dotplot(em2_new, showCategory=10, split=".sign")+ggtitle(paste0(cell,"_","gse_msigdbr_dotplot"))+facet_grid(.~.sign)
print(p)
dev.off()

}

set.seed(12345)

em2 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              pvalueCutoff = 0.3,
              verbose      = FALSE, seed=TRUE,eps=0, nPermSimple = 100000
              )
if(dim(as.data.frame(em2))[1]!=0){
em2=simplify(em2)
em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","GO_BP_gse"),quote=F,sep="\t")
#em2_new=read.table(paste0(cell,"_","GO_BP_gse"),header=T,sep="\t")

pdf(paste0(cell,"_","gse_GO_BP_dotplot","_top10.pdf"),width=8, height=7)
p=dotplot(em2_new, showCategory=10, split=".sign")+ggtitle(paste0(cell,"_","gse_GO_BP_dotplot"))+facet_grid(.~.sign)
print(p)
dev.off()
}

set.seed(12345)


em2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 10,
               pvalueCutoff = 0.1,
               verbose      = FALSE, seed=TRUE, eps=0, nPermSimple = 100000)
if(dim(as.data.frame(em2))[1]!=0){
em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","KEGG_gse"),quote=F,sep="\t")
#em2_new=read.table(paste0(cell,"_","KEGG_gse"),header=T,sep="\t")
pdf(paste0(cell,"_","gse_KEGG_dotplot","_top10.pdf"),width=8, height=7)
p=dotplot(em2_new, showCategory=10, split=".sign")+ggtitle(paste0(cell,"_","gse_KEGG_dotplot")) +facet_grid(.~.sign)
print(p)
dev.off()
}
#}
