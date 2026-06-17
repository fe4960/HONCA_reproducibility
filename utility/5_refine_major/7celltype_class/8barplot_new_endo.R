library(ggplot2)
dir="/dfs3b/ruic20_lab/junw42/"
cell="Endothelial_cell"
name="Endothelial_cell_subclass_sb_seurat_rmRPE"

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue.csv"),header=T,sep=",")


#level=c("Fibro_RPEchoroid_BMP5+", "Fibro_epipial", "Fibro_arachnoid_trabeculae_HMCN1+", "Fibro_arachnoid_barrier_STXBP6+", "Fibro_intima_pia", "Fibro_arachnoid_CLDN11+", "Fibro_dura_KCNMA1+", "Fibro_dura_boarder_SLC47A1+","Fibro_dura_NOX4+", "Fibro_arachnoid_TRPM3+",  "Fibro_x")


p=ggplot(data=data,aes(x=subclass,y=nCount_RNA,fill=tissue))+geom_bar(stat = "identity",position="fill")+labs(x = "Subclass", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("gold", "blue"))+xlab("")


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_tissue_df.pdf"))
print(p)
dev.off()

data=read.table(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source.csv"),header=T,sep=",")


p=ggplot(data=data,aes(x=subclass, y=nCount_RNA,fill=source))+geom_bar(stat = "identity",position="fill")+labs(x = "Data source", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))


pdf(paste0(dir,"HCA_ON/data/5_refine_major/scvi/",cell,"/clean/",name,"_source_df.pdf"))
print(p)
dev.off()


