args <- commandArgs(trailingOnly = TRUE)
library(qvalue)

cells=c("Astro_ONH","Astro_ONONH","Astro_retina","Astro_ON")
for(cell in cells){
dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Astrocyte_subclass_DEG/",cell)
setwd(dirIn)
fitmm=readRDS(paste0(cell,"_fitmm.rds"))
res=topTable( fitmm, coef='raceBlack', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_raceBlack"),sep="\t",quote=F)

#######res=topTable( fitmm, coef='raceAsian', number=Inf )
#######res$qval=qvalue(res$P.Value)$qvalues
######write.table(res,file=paste0(cell,"_DEG_raceAsian"),sep="\t",quote=F)

res=topTable( fitmm, coef='raceHispanic', number=Inf )
res$qval=qvalue(res$P.Value)$qvalues
write.table(res,file=paste0(cell,"_DEG_raceHispanic"),sep="\t",quote=F)
}
