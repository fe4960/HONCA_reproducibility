#peak=readRDS("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/peakset_clusters_major_full.rds")
peak=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/peakset_clusters_major_full.rds")
peak_anno=data.frame(peak)
peak_anno$peak=paste0(peak_anno$seqnames,":",peak_anno$start,"-",peak_anno$end)
#ca=read.table("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/cA_final_major_full_max",header=T,sep="\t")
ca=read.table("//dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/co_accessibility/ON_ATAC.txt.gz",header=T,sep="\t")

ca$peak1=sub("_",":",ca$queryPeak)  #paste0(ca$seqnames.1,":",ca$start,"-",ca$end)
ca$peak1=sub("_","-",ca$peak1) 

ca$peak2=sub("_",":",ca$subjectPeak)  #paste0(ca$seqnames.1,":",ca$start,"-",ca$end)
ca$peak2=sub("_","-",ca$peak2) 


#ca$peak2=#paste0(ca$seqnames.2,":",ca$start.1,"-",ca$end.1)
n=match(ca$peak1,peak_anno$peak,nomatch=0)
ca$peak1_nearestGene=peak_anno[n,]$nearestGene
ca$peak1_peakType=peak_anno[n,]$peakType
ca$peak1_distToGeneStart=peak_anno[n,]$distToGeneStart
ca$peak1_distToTSS=peak_anno[n,]$distToTSS
m=match(ca$peak2,peak_anno$peak,nomatch=0)


ca$peak2_nearestGene=peak_anno[m,]$nearestGene
ca$peak2_peakType=peak_anno[m,]$peakType
ca$peak2_distToGeneStart=peak_anno[m,]$distToGeneStart
ca$peak2_distToTSS=peak_anno[m,]$distToTSS

write.table(ca,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/cA_final_major_full_max_peakAnno",sep="\t",quote=F,row.names=F)

write.table(peak_anno,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/peakset_clusters_major_full_peak_anno",sep="\t",quote=F,row.names=F)

ca_enh=ca[ca$peak1_peakType=="Promoter"|ca$peak2_peakType=="Promoter",]
write.table(ca_enh,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/cA_final_major_full_max_peakAnno_LCRE",sep="\t",quote=F,row.names=F)

g2p=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/CRE/ON_ATAC_RNA_10000.txt.gz",header=T,sep="\t")

g2p$peak=sub("_", ":", g2p$peakName)
g2p$peak=sub("_", "-", g2p$peak)

#paste0(g2p$seqnames,":",g2p$start,"-",g2p$end)
LCRE=g2p[,c("geneName","peak")]
colnames(LCRE)=c("name","peak")
#LCRE=g2p[,c("name","peak")]
LCRE_ca1=ca_enh[ca$peak1_peakType=="Promoter",c("peak1_nearestGene","peak2")]
LCRE_ca2=ca_enh[ca$peak2_peakType=="Promoter",c("peak2_nearestGene","peak1")]
colnames(LCRE_ca1)=c("name","peak")
colnames(LCRE_ca2)=c("name","peak")

LCRE_all=rbind(LCRE,LCRE_ca1,LCRE_ca2)
LCRE_all_u=unique(LCRE_all)
# dim(LCRE_all)
#[1] 490567      2
# dim(LCRE_all_u)
#[1] 162492      2
write.table(LCRE_all_u, "/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/g2p_cA_final_major_full_max_LCRE_uniq",sep="\t",quote=F,row.names=F)
