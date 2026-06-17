library(ggplot2)
prop=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene_prop_age",sep="\t",header=T)

ggplot(data=prop,aes(x=factor(age),y=putative_sen))+geom_boxplot()+facet_wrap(~celltype)+ylim(0,0.25)
