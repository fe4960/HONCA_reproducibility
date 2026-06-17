library(motifStack)
library(motifbreakR)
data(motifbreakR_motif)

motif=new.env()
#ENCODE-motif    SOX9_known_ref8

#motif$NFE2L2=query(motifbreakR_motif,"MBD2_si")$'Hsapiens-HOCOMOCO-MBD2_si'
motif$NFE2L2=query(motifbreakR_motif,"SOX9_known_ref8")$'Hsapiens-ENCODE-motifs-SOX9_known_ref8'

#ref=read.table(file="sc_human_retina/scripts/ATAC_data/ref_MBD2",header=T,row.names=1)
ref=read.table(file="HCA_ON/data/13_GWAS/motif/ref_SOX9_rs7636574",header=T,row.names=1)

colnames(ref)=seq(1:32)

#alt=read.table(file="sc_human_retina/scripts/ATAC_data/alt_MBD2",header=T,row.names=1)
alt=read.table(file="HCA_ON/data/13_GWAS/motif/alt_SOX9_rs7636574",header=T,row.names=1)
colnames(alt)=seq(1:32)

motif$ref=data.matrix(ref)
motif$alt=data.matrix(alt)




motifs1 <- as.list(motif)

motifs4 <- mapply(motifs1, names(motifs1), FUN=function(.ele, .name){
  new("pfm",mat=.ele, name=.name)}, SIMPLIFY = FALSE)
hc <- clusterMotifs(motifs4)
library(ade4)
phylog <- ade4::hclust2phylog(hc)

pfmsAligned <- DNAmotifAlignment(motifs4)
## plot motifs
library(RColorBrewer)
color <- brewer.pal(5, "Set3")
#pdf("/storage/chenlab/Users/junwang/human_meta/data/finemap/MBD2.pdf")
pdf("HCA_ON/data/13_GWAS/motif/SOX9_rs7636574.pdf")

motifPiles(phylog=phylog, pfms=pfmsAligned,  #pfms=motifs4,  
            col.tree=rep(color, each=5),
            col.leaves=rep(rev(color), each=5),
#            col.pfms2=gpCol, 
#            r.anno=rep(0.02, length(dl)), 
#            col.anno=dl,
            motifScale="logarithmic",
            plotIndex=TRUE,
            groupDistance=10)
dev.off()
