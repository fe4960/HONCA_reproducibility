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


data1=rownames(data[(data$adj.P.Val<0.05)&(!is.na(data$adj.P.Val)),])


all_data=data1

genename=unique(all_data)


all_data_up=unique(rownames(data[(data$logFC>0)&(data$adj.P.Val<0.05)&(!is.na(data$adj.P.Val)),]))
all_data_down=unique(rownames(data[(data$logFC<0)&(data$adj.P.Val<0.05)&(!is.na(data$adj.P.Val)),]))


label=c("up","down")



library(data.table)
library('biomaRt')
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)

exp_gene_all=bitr(datalist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
set=list(all_data_up,all_data_down)

for(i in 1:length(set)){
genename=set[[i]]
label1=label[i]
eg_hm=bitr(genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gene=eg_hm$ENTREZID

library(msigdbr)
cate=c("C5")
for(c in cate){
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
m_t2g <- msigdbr(species = "Homo sapiens", category = c , subcategory="GO:BP") %>%
  dplyr::select(gs_name, entrez_gene)
em <- enricher(gene, TERM2GENE=m_t2g,minGSSize=10, qvalueCutoff =0.05, universe= exp_gene_all$ENTREZID, pvalueCutoff  = 0.05 , pAdjustMethod = "BH")
if(dim(as.data.frame(em))[1]!=0){
#####write.table(data.frame(em),file=paste0(cell,"_msigdbr_",c,"_enricher_minGSSize10","_",label1),quote=F,sep="\t")

em_new <- setReadable(em, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
write.table(data.frame(em_new),file=paste0(cell,"_msigdbr_",c,"_enricher","_",label1),quote=F,sep="\t")


pdf(paste0("region_",cell,"_","msigdbr_",c,"_enricher_minGSSize10_dotPlot_Cate20","_",label1,".pdf"),width=10,height=20)
print(dotplot(em, showCategory=20))
dev.off()


}
}



library(GOSemSim)
ego <- enrichGO(gene  = gene,
        universe      = exp_gene_all$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize = 10,
        readable      = TRUE)

if(dim(as.data.frame(ego))[1]!=0){

ego=simplify(ego)


    pdf(paste0("region_",cell,"_","GO_MF_dotplot","_",label1,".pdf"),width=10,height=20)

    print(dotplot(ego, showCategory = 20 ))

    dev.off()
    em_new <- setReadable(ego, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
    write.table(data.frame(em_new),file=paste0("enrichGO_MF_readable_",cell,"_",label1),quote=F,sep="\t")
}

ego <- enrichGO(gene  = gene,
        universe      = exp_gene_all$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize = 10,
        readable      = TRUE)
if(dim(as.data.frame(ego))[1]!=0){
ego=simplify(ego)

    pdf(paste0("region_",cell,"_","GO_BP_dotplot","_",label1,".pdf"),width=10,height=20)

	print(dotplot(ego, showCategory = 20 ))
    dev.off()
    em_new <- setReadable(ego, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
    write.table(data.frame(em_new),file=paste0("enrichGO_BP_readable_",cell,"_",label1),quote=F,sep="\t")
}

ego <- enrichKEGG(gene  = gene,
        universe      = exp_gene_all$ENTREZID,
        organism         = "hsa",
	keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        minGSSize = 10)
if(dim(as.data.frame(ego))[1]!=0){

    pdf(paste0("region_",cell,"_","KEGG_BP_dotplot","_",label1,".pdf"),width=10,height=20)

	print(dotplot(ego, showCategory = 20 ))
    dev.off()
    em_new <- setReadable(ego, OrgDb ="org.Hs.eg.db", keyType="ENTREZID")
    write.table(data.frame(em_new),file=paste0("enrichKEGG_readable_",cell,"_",label1),quote=F,sep="\t")
}
}


