library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
i=args[1]
fn=args[2]
cell=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/8sample_celltype")
#for(i in cell$V1){
file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/",i,"/clean/",fn,"_sample_number.txt")

#file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/",i,"/clean/",i,"_subclass_sample_number.txt")
data=read.table(file,header=T)
pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/",i,"/clean/",fn,"_sample_number.pdf"))

#pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/",i,"/clean/",i,"_subclass_sample_number.pdf"))
my_palette <- c("lightblue", "orange", "lightgreen","grey","blue", "red", "purple1","yellow","darkgreen","pink","darkblue","thistle","green","lightyellow")
colors_to_use <- rep(my_palette, ceiling(length(unique(data$sampleid)) / length(my_palette)))
p=ggplot(data,aes(x=subclass,y=nCount_RNA,fill=sampleid))+geom_bar(position="fill",stat="identity")+coord_flip()+xlab("")+ylab("")+theme(text = element_text(size=20),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"),legend.position="none")+scale_fill_manual(values = colors_to_use)+xlab("Sample distribution")
		#scale_color_brewer(palette = "Set1",direction=-1)
print(p)
dev.off()
#}
