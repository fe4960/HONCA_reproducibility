library(ggplot2)
cells=c("all","Oligodendrocyte","Astrocyte","Fibroblast","Endothelial_cell")
for(cell in cells){
data=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",cell,"_subclass_DEG/",cell,"_DEG_p_fc_padj0.05"),sep="\t")
data=data[(data$V2!="AC")&(data$V2!="Rod")&(data$V2!="BC"),]
data$label="Up"
data[data$V3<0,"label"]="Down"
p=ggplot(data=data,aes(x=V2,y=V3,fill=label))+geom_boxplot()+theme(text = element_text(size=25),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"), axis.text.x = element_text(angle = 90, hjust = 1))+ylab("Gene expression Log2FC")+xlab("")+scale_fill_manual(values= c("Up" = "gold", "Down" = "purple"))+geom_hline(yintercept = 0, linetype = "dashed")+ggtitle(cell)
pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",cell,"_subclass_DEG/",cell,"_DEG_p_fc_all.pdf"),height=8)
print(p)
dev.off()
}

