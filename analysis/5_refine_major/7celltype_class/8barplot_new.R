library(ggplot2)
dir="/dfs3b/ruic20_lab/junw42/"
#cell="Astrocyte"
#name="Astrocyte_subclass_new5_clean"
cell="Fibroblast"
name="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location"

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue.csv"),header=T,sep=",")

#level=c("Astro_ON1_WNK2+TSHZ2+","Astro_ON2_ACTN1+SERPINA3+","Astro_ON3_NR4A3+FOS+","Astro_ON4_AFF3+DMGDH+","Astro_ON5_DPP10+","Astro_ONHON1_SLC4A11+MARCH1+","Astro_ONHON2_CST3+APOE+","Astro_ONH1_GABBR2+SV2B+","Astro_ONH2_PAX5+GABBR2+","Astro_retina1_PAX5+ME1+","Astro_retina2_NLGN1+")

level=c("Fibro_RPEchoroid_BMP5+", "Fibro_epipial", "Fibro_arachnoid_trabeculae_HMCN1+", "Fibro_arachnoid_barrier_STXBP6+", "Fibro_intima_pia", "Fibro_arachnoid_CLDN11+", "Fibro_dura_KCNMA1+", "Fibro_dura_boarder_SLC47A1+","Fibro_dura_NOX4+", "Fibro_arachnoid_TRPM3+",  "Fibro_x")

#{"0":"Fibro_intima_pia" ,"13":"Fibro_intima_pia" ,"8": "Fibro_epipial","6":"Fibro_RPEchoroid_BMP5+","7":"Fibro_dura_NOX4+","3":"Fibro_dura_KCNMA1+","1":"Fibro_dura_boarder_SLC47A1+","5":"Fibro_arachnoid_barrier_STXBP6+","2":"Fibro_arachnoid_CLDN11+","4": "Fibro_arachnoid_TRPM3+","10":"Fibro_x", "9": "Fibro_arachnoid_trabeculae_HMCN1+"}

p=ggplot(data=data,aes(x=factor(subclass,levels=level),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")

#p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue_df.pdf"))
print(p)
dev.off()

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source.csv"),header=T,sep=",")


p=ggplot(data=data,aes(x=factor(subclass, levels=level),y=nCount_RNA,fill=source))+geom_bar(stat = "identity",position="fill")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source_df.pdf"))
print(p)
dev.off()

cell="Mural_cell"
name="Mural_cell_subclass_rmRPE"

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue.csv"),header=T,sep=",")

#level=c("Astro_ON1_WNK2+TSHZ2+","Astro_ON2_ACTN1+SERPINA3+","Astro_ON3_NR4A3+FOS+","Astro_ON4_AFF3+DMGDH+","Astro_ON5_DPP10+","Astro_ONHON1_SLC4A11+MARCH1+","Astro_ONHON2_CST3+APOE+","Astro_ONH1_GABBR2+SV2B+","Astro_ONH2_PAX5+GABBR2+","Astro_retina1_PAX5+ME1+","Astro_retina2_NLGN1+")

level=c("Fibro_RPEchoroid_BMP5+", "Fibro_epipial", "Fibro_arachnoid_trabeculae_HMCN1+", "Fibro_arachnoid_barrier_STXBP6+", "Fibro_intima_pia", "Fibro_arachnoid_CLDN11+", "Fibro_dura_KCNMA1+", "Fibro_dura_boarder_SLC47A1+","Fibro_dura_NOX4+", "Fibro_arachnoid_TRPM3+",  "Fibro_x")

p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")



pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue_df.pdf"))
print(p)
dev.off()

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source.csv"),header=T,sep=",")


p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=source))+geom_bar(stat = "identity",position="fill")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source_df.pdf"))
print(p)
dev.off()


cell="Macrophage"
name="Macrophage_subclass_sb_clean_rmRPE"

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue.csv"),header=T,sep=",")


level=c("Fibro_RPEchoroid_BMP5+", "Fibro_epipial", "Fibro_arachnoid_trabeculae_HMCN1+", "Fibro_arachnoid_barrier_STXBP6+", "Fibro_intima_pia", "Fibro_arachnoid_CLDN11+", "Fibro_dura_KCNMA1+", "Fibro_dura_boarder_SLC47A1+","Fibro_dura_NOX4+", "Fibro_arachnoid_TRPM3+",  "Fibro_x")

p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")



pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue_df.pdf"))
print(p)
dev.off()

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source.csv"),header=T,sep=",")


p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=source))+geom_bar(stat = "identity",position="fill")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source_df.pdf"))
print(p)
dev.off()


cell="Microglia"
name="Microglia_subclass_sb_clean"

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue.csv"),header=T,sep=",")


level=c("Fibro_RPEchoroid_BMP5+", "Fibro_epipial", "Fibro_arachnoid_trabeculae_HMCN1+", "Fibro_arachnoid_barrier_STXBP6+", "Fibro_intima_pia", "Fibro_arachnoid_CLDN11+", "Fibro_dura_KCNMA1+", "Fibro_dura_boarder_SLC47A1+","Fibro_dura_NOX4+", "Fibro_arachnoid_TRPM3+",  "Fibro_x")

p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")



pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue_df.pdf"))
print(p)
dev.off()

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source.csv"),header=T,sep=",")


p=ggplot(data=data,aes(x=factor(subclass),y=nCount_RNA,fill=source))+geom_bar(stat = "identity",position="fill")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source_df.pdf"))
print(p)
dev.off()
