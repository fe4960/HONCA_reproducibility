
library(ggalluvial)
library(dplyr)
library(ggplot2)


plot_func=function(file){
vac=read.table(file=file, sep=",", header=T)

colnames(vac)=c("PMID37566633","This_study", "count")

#sum1=sum(as.numeric(vac[,"count"]))

vac$count=as.numeric(vac$count)
vac$freq=vac$count/sum(vac$count)


p=ggplot(data = vac,
       aes(axis1 = PMID37566633, axis2 = This_study, y = freq)) +
  geom_alluvium(aes(fill = This_study)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("PMID37566633", "This_study"),
                   expand = c(0.15, 0.05)) +
  theme_void() 


pdf1=paste0(indir,"/", bname,"_sanes_compare.pdf") #, height=12)

pdf(pdf1)
print(p)
dev.off()

}

args <- commandArgs(trailingOnly = TRUE)
dir="/dfs3b/ruic20_lab/junw42/"
bname=args[1]
celltype=args[2]
indir=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/",celltype)
#outdir=f"{indir}/clean"


#file=paste0(indir,"/clean/", bname,"_sanes_compare.csv")
file=paste0(indir,"/clean/", bname,"_sanes_compare_edit.csv")

plot_func(file) 

