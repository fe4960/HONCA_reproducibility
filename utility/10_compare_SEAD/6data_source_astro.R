library(ggplot2)
dir="/dfs3b/ruic20_lab/junw42/"
data=read.table(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.3_leiden_disease"),header=T,sep=",")
pdf(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.3_leiden_disease.pdf"))
ggplot(data,aes(x=disease,y=nCount_RNA,fill=factor(leiden1)))+geom_bar(position="fill",stat="identity")+coord_flip()+scale_fill_brewer(palette="Set2")+theme_minimal()+theme(text=element_text(size=20))+labs(fill="leiden\ncluster")
dev.off()


data=read.table(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.3_leiden_disease_ONONH_PFC"),header=T,sep=",")
pdf(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.3_leiden_disease_ONONH_PFC.pdf"))
ggplot(data,aes(x=disease,y=nCount_RNA,fill=factor(leiden1)))+geom_bar(position="fill",stat="identity")+coord_flip()+scale_fill_brewer(palette="Set2")+theme_minimal()+theme(text=element_text(size=20))+labs(fill="leiden\ncluster")
dev.off()


data=read.table(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.4_leiden_disease"),header=T,sep=",")
pdf(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.4_leiden_disease.pdf"))
ggplot(data,aes(x=disease,y=nCount_RNA,fill=factor(leiden1)))+geom_bar(position="fill",stat="identity")+coord_flip()+scale_fill_brewer(palette="Set2")+theme_minimal()+theme(text=element_text(size=20))+labs(fill="leiden\ncluster")
dev.off()


data=read.table(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.4_leiden_disease_ONONH_PFC"),header=T,sep=",")
pdf(paste0(dir,"HCA_ON/data/scvi/RNA/raw/Astrocyte_GSE267301_merge/rs_0.4_leiden_disease_ONONH_PFC.pdf"))
ggplot(data,aes(x=disease,y=nCount_RNA,fill=factor(leiden1)))+geom_bar(position="fill",stat="identity")+coord_flip()+scale_fill_brewer(palette="Set2")+theme_minimal()+theme(text=element_text(size=20))+labs(fill="leiden\ncluster")
dev.off()

