library(UpSetR)
cells=c("major") #Fibroblast","Endothelial_cell","Astrocyte","Oligodendrocyte","all")
for(cell in cells){
file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",cell,"_majorclass_DEG/",cell,"_upset_list_age_cpm1_all_ct_nomashr_tab")
data=read.table(file,header=F,sep="\t")

#input=c(t(data))
#strsplit("=",input)
input=data$V2
names(input)=data$V1
pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",cell,"_majorclass_DEG/",cell,"_upset_list_age_cpm1_all_ct_nomashr.pdf"),width=10)

p=upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = cell,
      sets.bar.color="blue"
      )
print(p)

dev.off()
}

