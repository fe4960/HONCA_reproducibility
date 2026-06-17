library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)

ct=args[1]
meta=args[2]
fn=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_",ct)
#wd==paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_",ct)
data=read.table(fn,header=T,row.names=1,sep=",")
file_info=read.table(meta,header=T,sep=",")

fListNames=colnames(data)
fListNames=gsub("X","",fListNames)
fListNames=gsub("\\.","\\-",fListNames)

m=match(fListNames,file_info$donor,nomatch=0)

metadata=file_info[m,]
metadata$age=as.numeric(metadata$age_year)

#data1=t(data2)
#DT1=as.matrix(scale(data1,center=T))
DT1=as.matrix(scale(data,center=F))
rownames(DT1)=rownames(data)
colnames(DT1)=colnames(data)

my_sample_col <- data.frame(age = metadata$age,pt=metadata$donor)

rownames(my_sample_col) <- colnames(DT1)

my_sample_col=my_sample_col[order(my_sample_col$age),]

#n=match(rownames(DT1),meta$donor,nomatch=0)

#my_sample_row = data.frame(celltype=rownames(DT1))

#my_sample_row$age_range="60younger"

#my_sample_row[my_sample_row$age>=60,]$age_range="60orOlder"

rownames(my_sample_row) = rownames(DT1)


#my_sample_row=my_sample_row[order(my_sample_row$age),]

n=match(names(my_sample_col), colnames(DT1),nomatch=0)

DT1=DT1[,n]

#ha=HeatmapAnnotation(Celltype=my_sample_col[colnames(DT1),]$celltype,Exp_Dir=my_sample_col[colnames(DT1),]$dir, col=list(Celltype=c("OLIGO"="steelblue","OPC"="grey", "Macro"="darkgreen", "Micro"="yellow", "MG"="lightgreen","Astro"="thistle"), Exp_Dir=c("Down"="skyblue","Up"="salmon")),annotation_name_side = "right",annotation_name_gp= gpar(fontsize = 15))
#ha=HeatmapAnnotation(age=my_sample_col$age,annotation_name_side = "right",annotation_name_gp= gpar(fontsize = 15))
#va=rowAnnotation(Age_range=my_sample_row[rownames(DT1),]$age_range, Age=my_sample_row[rownames(DT1),]$age ,  col=list(Age_range=c("60younger"="gold","60orOlder"="purple"))) 

library(ggplot2)
#setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_noPlot/")
outfile=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final/snRNA_donor_num_",ct,"_prop_heatmap.pdf")
#pdf(paste0("8heatmap_",m,"_glia_only_nolog_rmMG_age.pdf"),width=12,height=6)
pdf(outfile, width=12,height=6)
#ht=Heatmap(DT1,row_title = paste0(ct),  clustering_distance_rows = "euclidean", show_row_names = T, name = paste0(ct),  top_annotation = ha, cluster_columns = F, cluster_rows = T) #,  left_annotation=va)
#draw(ht)

#ggplot(data=DT1,aes(x=))
dev.off()

