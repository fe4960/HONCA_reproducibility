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
ct=args[1]
cell=args[2]
#seq=args[2]
#dem=args[3]
j=as.numeric(args[3])
dem=args[4]
#dirIn=paste0("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/",cell,"_",seq,"_",dem)
#dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",ct,"_subclass_DEG/",cell,"/",cell,"_DEG_GO_",dem)
dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",ct,"_majorclass_DEG/",cell,"/",cell,"_DEG_GO_",dem)
dir.create(dirIn,recursive = TRUE)
setwd(dirIn)
dirIn <- getwd()
all_data=NULL

set.seed(1234567)
#file=paste0("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/",dem,"_DEG_res_cpm_",dem,"_pm_ct_ancFix")
#file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",ct,"_subclass_DEG/",dem,"_DEG_res_cpm_",dem,"_pm_ct_ancFix")
#####file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",ct,"_subclass_DEG/",dem,"_DEG_res_cpm_",dem)
fl=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/", ct, "_majorclass_DEG/",cell,"/", cell, "_DEG_", dem)

data=read.table(fl,header=T)

#fl="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas_sample/cpm_1_gene_ct.txt"

#fl=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",ct,"_subclass_DEG/",dem,"_DEG_res_cpm_",dem,"_beta_ancFix")

genel=read.table(fl,header=T)



datalist=rownames(genel) #rownames(data)

engene_entrez= bitr(datalist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

n=match(engene_entrez$SYMBOL,datalist,nomatch=0)

geneList=data[n,j]


m=match(datalist[n],engene_entrez$SYMBOL,nomatch=0)


names(geneList)=engene_entrez[m,]$ENTREZID


geneList=sort(geneList, decreasing = TRUE)

unique_geneList <- geneList[!duplicated(names(geneList))]
geneList=sort(unique_geneList, decreasing = TRUE)

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
write.table(data.frame(em2_new),file=paste0(cell,"_","msigdbr_gse","_",dem),quote=F,sep="\t")
pdf(paste0(cell,"_","gsea_msigdbr_dotplot_",dem,".pdf"),width=18, height=20)
p=dotplot(em2_new, showCategory=20, split=".sign")+ggtitle(paste0(cell,"_","gse_msigdbr_dotplot_",dem))+facet_grid(.~.sign)
print(p)
dev.off()

}

set.seed(12345)

em2 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              pvalueCutoff = 0.1,
              verbose      = FALSE, seed=TRUE,eps=0, nPermSimple = 100000
              )
if(dim(as.data.frame(em2))[1]!=0){
em2=simplify(em2)
em2_new <- setReadable(em2, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em2_new),file=paste0(cell,"_","GO_BP_gse","_",dem),quote=F,sep="\t")
pdf(paste0(cell,"_","gse_GO_BP_dotplot_",dem,".pdf"),width=18, height=20)
p=dotplot(em2_new, showCategory=20, split=".sign")+ggtitle(paste0(cell,"_","gse_GO_BP_dotplot_",dem))+facet_grid(.~.sign)
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
write.table(data.frame(em2_new),file=paste0(cell,"_","KEGG_gse","_",dem),quote=F,sep="\t")
pdf(paste0(cell,"_","gse_KEGG_dotplot_",dem,".pdf"),width=18, height=20)
p=dotplot(em2_new, showCategory=20, split=".sign")+ggtitle(paste0(cell,"_","gse_KEGG_dotplot_",dem)) +facet_grid(.~.sign)
print(p)
dev.off()

}

