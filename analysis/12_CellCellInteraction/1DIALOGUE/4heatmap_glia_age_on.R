library("ComplexHeatmap")
indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass"
oligo=read.table(paste0(indir,"/exp_Oligodendrocyte_logNorm_all"),header=T)
opc=read.table(paste0(indir,"/exp_Oligodendrocyte_precursor_cell_logNorm_all"),header=T)
#macro=read.table(paste0(indir,"/exp_Macrophage_logNorm_all"),header=T)
astro=read.table(paste0(indir,"/exp_Astrocyte_logNorm_all"),header=T)
micro=read.table(paste0(indir,"/exp_Microglia_logNorm_all"),header=T)
#setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot/")
#setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_2k/")
#setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_2k_noAge/")
#setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_2k_Age/")
#setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_3k_Age/")
setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_3k_AgeRange/")


#mg=read.table(paste0(indir,"/exp_MG_logNorm_all"),header=T)

########oligo=read.table(paste0(indir,"/exp_Oligodendrocyte_Norm_all"),header=T)
#########opc=read.table(paste0(indir,"/exp_Oligodendrocyte_precursor_cell_Norm_all"),header=T)
#macro=read.table(paste0(indir,"/exp_Macrophage_logNorm_all"),header=T)
########astro=read.table(paste0(indir,"/exp_Astrocyte_Norm_all"),header=T)
########micro=read.table(paste0(indir,"/exp_Microglia_Norm_all"),header=T)
#########mg=read.table(paste0(indir,"/exp_MG_Norm_all"),header=T)


#indir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass/Astrocyte_subclass"

#astro_ret=read.table(paste0(indir1,"/exp_Astro_ret_logNorm_all"),header=T)

#astro_on=read.table(paste0(indir1,"/exp_Astro_ON_logNorm_all"),header=T)

common_cols=Reduce(intersect, list(colnames(oligo), colnames(opc),  colnames(astro)  , colnames(micro)))

#####common_cols=Reduce(intersect, list(colnames(oligo), colnames(opc),  colnames(astro)  , colnames(micro), colnames(mg)))
meta=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta",sep=",",header=T)
#data=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE//glia_only_on_noPlot/DLG.full.output_ONONH_glia_only.rds")
data=readRDS("DLG.full.output_ONONH_glia_only.rds")
#mcp=c("MCP2","MCP4","MCP5")
mcp=c("MCP1","MCP2","MCP3")
for(m in mcp){

oligo_up=oligo[data$MCPs[[m]]$oligo.up,common_cols]
oligo_down=oligo[data$MCPs[[m]]$oligo.down,common_cols]

opc_up=opc[data$MCPs[[m]]$opc.up,common_cols]
opc_down=opc[data$MCPs[[m]]$opc.down,common_cols]

######oligo_up=oligo[data$MCPs.full[[m]]$oligo.up,common_cols]
#######oligo_down=oligo[data$MCPs.full[[m]]$oligo.down,common_cols]

######opc_up=opc[data$MCPs.full[[m]]$opc.up,common_cols]
#######opc_down=opc[data$MCPs.full[[m]]$opc.down,common_cols]

#macro_up=oligo[data$MCPs[[m]]$macro.up,common_cols]
#macro_down=oligo[data$MCPs[[m]]$macro.down,common_cols]

######micro_up=micro[data$MCPs.full[[m]]$micro.up,common_cols]
#######micro_down=micro[data$MCPs.full[[m]]$micro.down,common_cols]

micro_up=micro[data$MCPs[[m]]$micro.up,common_cols]
micro_down=micro[data$MCPs[[m]]$micro.down,common_cols]


#######mg_up=mg[data$MCPs.full[[m]]$MG.up,common_cols]
#######mg_down=mg[data$MCPs.full[[m]]$MG.down,common_cols]

#####astro_up=astro[data$MCPs.full[[m]]$astro.up,common_cols]
#####astro_down=astro[data$MCPs.full[[m]]$astro.down,common_cols]


astro_up=astro[data$MCPs[[m]]$astro.up,common_cols]
astro_down=astro[data$MCPs[[m]]$astro.down,common_cols]

#astroret_up=oligo[data$MCPs[[m]]$astroret.up,common_cols]
#astroret_down=oligo[data$MCPs[[m]]$astroret.down,common_cols]

#####data2=rbind(oligo_up,opc_up, astro_up, micro_up, mg_up,  oligo_down, opc_down, astro_down,  micro_down,  mg_down)
data2=rbind(oligo_up,opc_up, astro_up, micro_up,  oligo_down, opc_down, astro_down,  micro_down)

data1=t(data2)
#DT1=as.matrix(t(scale(t(data1),center=T)))
DT1=as.matrix(scale(data1,center=T))
#######my_sample_col <- data.frame(celltype =c(rep("OLIGO",dim(oligo_up)[1]), rep("OPC",dim(opc_up)[1]), rep("Astro",dim(astro_up)[1]), rep("Micro",dim(micro_up)[1]),  rep("MG",dim(mg_up)[1]),  rep("OLIGO",dim(oligo_down)[1]),  rep("OPC",dim(opc_down)[1]), rep("Astro",dim(astro_down)[1]), rep("Micro",dim(micro_down)[1]), rep("MG",dim(mg_down)[1]) ), dir=c(rep("Up",dim(oligo_up)[1]), rep("Up",dim(opc_up)[1]),  rep("Up",dim(astro_up)[1]),   rep("Up",dim(micro_up)[1]), rep("Up",dim(mg_up)[1]),  rep("Down",dim(oligo_down)[1]), rep("Down",dim(opc_down)[1]),  rep("Down",dim(astro_down)[1]),rep("Down",dim(micro_down)[1]), rep("Down",dim(mg_down)[1]) ))

my_sample_col <- data.frame(celltype =c(rep("OLIGO",dim(oligo_up)[1]), rep("OPC",dim(opc_up)[1]), rep("Astro",dim(astro_up)[1]), rep("Micro",dim(micro_up)[1]),    rep("OLIGO",dim(oligo_down)[1]),  rep("OPC",dim(opc_down)[1]), rep("Astro",dim(astro_down)[1]), rep("Micro",dim(micro_down)[1]) ), dir=c(rep("Up",dim(oligo_up)[1]), rep("Up",dim(opc_up)[1]),  rep("Up",dim(astro_up)[1]),   rep("Up",dim(micro_up)[1]),  rep("Down",dim(oligo_down)[1]), rep("Down",dim(opc_down)[1]),  rep("Down",dim(astro_down)[1]),rep("Down",dim(micro_down)[1]) ))



row.names(my_sample_col) <- colnames(DT1)

n=match(rownames(DT1),meta$donor,nomatch=0)

my_sample_row = data.frame(age=meta[n,]$age_year)

my_sample_row$age_range="60younger"

my_sample_row[my_sample_row$age>=60,]$age_range="60orOlder"

row.names(my_sample_row) = rownames(DT1)

my_sample_row=my_sample_row[order(my_sample_row$age),]

n=match(rownames(my_sample_row), rownames(DT1),nomatch=0)

DT1=DT1[n,]

ha=HeatmapAnnotation(Celltype=my_sample_col[colnames(DT1),]$celltype,Exp_Dir=my_sample_col[colnames(DT1),]$dir, col=list(Celltype=c("OLIGO"="steelblue","OPC"="lightgreen", "Macro"="darkgreen", "Micro"="yellow", "MG"="grey","Astro"="thistle"), Exp_Dir=c("Down"="skyblue","Up"="salmon")),annotation_name_side = "right",annotation_name_gp= gpar(fontsize = 15))

#va=rowAnnotation(Age_range=my_sample_row[rownames(DT1),]$age_range,  col=list(Age_range=c("60younger"="gold","60orOlder"="purple"))) #,name_side = "right",name_gp= gpar(fontsize = 15))

va=rowAnnotation(Age_range=my_sample_row[rownames(DT1),]$age_range, Age=my_sample_row[rownames(DT1),]$age ,  col=list(Age_range=c("60younger"="gold","60orOlder"="purple"))) #,name_side = "right",name_gp= gpar(fontsize = 15))


####pdf(paste0("8heatmap_",m,"_glia_only_log_full_rmMG.pdf"),width=12,height=6)
#pdf(paste0("8heatmap_",m,"_glia_only_log_rmMG.pdf"),width=12,height=6)
#pdf(paste0("8heatmap_",m,"_glia_only_nolog_rmMG.pdf"),width=12,height=6)
#pdf(paste0("8heatmap_",m,"_glia_only_on_log_age_spearman.pdf"),width=12,height=6)
#pdf(paste0("8heatmap_",m,"_glia_only_on_log_age_pearson.pdf"),width=12,height=6)
#pdf(paste0("8heatmap_",m,"_glia_only_on_log_noage_pearson.pdf"),width=12,height=6)
#pdf(paste0("8heatmap_",m,"_glia_only_on_log_noage_euclidean.pdf"),width=12,height=6)
pdf(paste0("8heatmap_",m,"_glia_only_on_log_noage_spearman.pdf"),width=12,height=6)


#clustering_distance_rows = "spearman" #"euclidean",
ht=Heatmap(DT1,row_title = paste0(m,"_glia"),  clustering_distance_rows = "spearman", clustering_method_rows="ward.D2", show_row_names = FALSE, name = paste0(m,"_glia"),  top_annotation = ha, cluster_columns = F, cluster_rows = T,  left_annotation=va)
draw(ht)
dev.off()
}


#n=match(c
#oligo1=oligo[,(colnames(oligo) %in% colnames(opc))&&(colnames(oligo) %in% colnames(macro))&& (colnames(oligo) %in% colnames(astro)) && (colnames(oligo) %in% colnames(micro)) && (colnames(oligo) %in% colnames(mg))  ]
