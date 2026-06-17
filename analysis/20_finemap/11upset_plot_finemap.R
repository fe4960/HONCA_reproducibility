library(UpSetR)
cells=c("major") #Fibroblast","Endothelial_cell","Astrocyte","Oligodendrocyte","all")
for(cell in cells){
#file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",cell,"_majorclass_DEG/",cell,"_upset_list_age_cpm1_all_ct_nomashr_tab")
file=paste0("/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/pip_var_all_merged_uniform_anno_all_ON_01")

#data=read.table(file,header=F,sep="\t")

#input=read.table(file,header=F,sep="\t")

input <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#input=c(t(data))
tmp=strsplit(input$V1, "&=" )
#########input=data$V2
#########names(input)=data$V1

#counts <- as.numeric(sapply(tmp, "[", 2))
#names(counts) <- sapply(tmp, "[", 1)

counts <- as.numeric(sapply(tmp, `[[`, 2))
names(counts) <- sapply(tmp, `[[`, 1)


#pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/",cell,"_majorclass_DEG/",cell,"_upset_list_age_cpm1_all_ct_nomashr.pdf"),width=10)
pdf(paste0("/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/pip_var_all_merged_uniform_anno_all_ON_01.pdf"),width=10)

p=upset(fromExpression(counts), 
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
      sets.x.label = "major",
      sets.bar.color="blue"
      )
print(p)

dev.off()
}

