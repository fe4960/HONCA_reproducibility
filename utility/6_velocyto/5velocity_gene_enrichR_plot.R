library("ggplot2")
library(dplyr)
library(tidyr)
library(stringr)
#library(enrichR)
#args <- commandArgs(trailingOnly = TRUE)


#dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023","GO_Biological_Process_2023","MSigDB_Hallmark_2020")
dbs <- c("MSigDB_Hallmark_2020","GO_Biological_Process_2023","GO_Molecular_Function_2023")

plot=function(dirIn,cell){
setwd(dirIn)

for(d in dbs){
data=read.table(paste0("enrichr_",d,".txt"),header=T,sep="\t")
head(data)
data=data %>% separate(Overlap,into=c("Gene_count","backgroud_count"),sep="/")
data=data %>% separate(Term,into=c("Term","other"),sep = "\\(|\\)")

#data=data[data$Adjusted.P.value<=0.05,]
data=data[1:20,]
data$Gene_count=as.numeric(data$Gene_count)
data_sorted = data %>% arrange(Gene_count)
p=ggplot(data_sorted,aes(x=Gene_count,y=factor(Term,levels=data_sorted$Term)))+geom_point(aes(size=as.numeric(Odds.Ratio), color=as.numeric(Adjusted.P.value)))+theme_minimal()+theme(text=element_text(size=25))+scale_size_continuous(range=c(2,15))+scale_color_continuous(low = "blue", high = "red")+xlab("Gene count")+labs(size="Odds.Ratio",color="Padj")+ylab("Term")+scale_y_discrete(labels = function(y) str_wrap(y, width = 40))

#ggplot(data=data,aes(x=V1,y=V2))+geom_point(aes(size=V3),color=data$sig)+theme_minimal()+theme(text=element_text(size=20),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(2,15))+scale_y_discrete(labels = function(y) str_wrap(y, width = 40))


pdf(paste0("enrichr_",d,"_",cell,".pdf"),width=12,height=15)
print(p)
dev.off()
}


}


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new_update_enrichR")
dir.create(dirIn,recursive = TRUE)

plot(dirIn,"Astrocyte")

dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/oligo_opc_hvg_2000_raw_new_enrichR")
dir.create(dirIn,recursive = TRUE)
plot(dirIn,"OPC_Oligo")


dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/veloc/Microglia_subclass_sb_hvg_2000_raw_enrichR")
dir.create(dirIn,recursive = TRUE)
plot(dirIn,"Microglia")

##########

dirIn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/veloc/Macrophage_subclass_sb_hvg_2000_raw_enrichR")
dir.create(dirIn,recursive = TRUE)
plot(dirIn,"Macrophage")

