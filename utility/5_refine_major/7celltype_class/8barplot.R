library(ggplot2)
dir="/dfs3b/ruic20_lab/junw42/"
#cell="Oligodendrocyte_precursor_cell"
#cell="Oligodendrocyte"
#name="Oligodendrocyte_subclass_seurat"
#name="Oligodendrocyte_precursor_cell_subclass_seurat"
cell="Astrocyte"
name="Astrocyte_subclass_new5_clean"
#cell="Fibroblast"
#name="Fibroblast_subclass_clean_new"
data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue.csv"),header=T,sep=",")

#data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_tissue_df.csv"),header=T,sep=",")
#p=ggplot(data=data,aes(x=subclass,y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 60, hjust = 1), text=element_text(size = 25))+scale_fill_manual(values = c("gold", "blue"))
####p=ggplot(data=data,aes(x=factor(subclass,levels=c("Astro_ON1_TSHZ2+","Astro_ON2_NR4A3+","Astro_ON3_AFF3+","Astro_ON4_TSHZ2+DPP10+","Astro_ONHON1_SLC4A11+","Astro_ONHON2_SLC4A11+APOE+","Astro_ONHON3_ANKUB1+","Astro_ONH1_GABBR2+","Astro_ONH2_PAX5+GABBR2+","Astro_retina1_PAX5+ME1+","Astro_retina2_NLGN1+")),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))
level=c("Astro_ON1_WNK2+TSHZ2+","Astro_ON2_ACTN1+SERPINA3+","Astro_ON3_NR4A3+FOS+","Astro_ON4_AFF3+DMGDH+","Astro_ON5_DPP10+","Astro_ONHON1_SLC4A11+MARCH1+","Astro_ONHON2_CST3+APOE+","Astro_ONH1_GABBR2+SV2B+","Astro_ONH2_PAX5+GABBR2+","Astro_retina1_PAX5+ME1+","Astro_retina2_NLGN1+")
p=ggplot(data=data,aes(x=factor(subclass,levels=level),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")

pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue_df.pdf"))
print(p)
dev.off()

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source.csv"),header=T,sep=",")

#data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_source_df.csv"),header=T,sep=",")

p=ggplot(data=data,aes(x=subclass,y=nCount_RNA,fill=source))+geom_bar(stat = "identity",position="fill")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))


#p=ggplot(data=data,aes(x=subclass,y=nCount_RNA,fill=source))+geom_bar(stat = "identity")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 60, hjust = 1), text=element_text(size = 25))+scale_fill_manual(values = c("yellow", "purple"))
#pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_source_df.pdf"))
pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source_df.pdf"))
print(p)
dev.off()
