library("ComplexHeatmap")
indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass"
oligo=read.table(paste0(indir,"/exp_Oligodendrocyte_logNorm_all"),header=T)
#oligo=oligo.T
opc=read.table(paste0(indir,"/exp_Oligodendrocyte_precursor_cell_logNorm_all"),header=T)
macro=read.table(paste0(indir,"/exp_Macrophage_logNorm_all"),header=T)
#astro=read.table(paste0(indir,"/exp_Astrocyte_logNorm_all"),header=T)
micro=read.table(paste0(indir,"/exp_Microglia_logNorm_all"),header=T)
mg=read.table(paste0(indir,"/exp_MG_logNorm_all"),header=T)

indir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass/Astrocyte_subclass"

astro_ret=read.table(paste0(indir1,"/exp_Astro_ret_logNorm_all"),header=T)

astro_on=read.table(paste0(indir1,"/exp_Astro_ON_logNorm_all"),header=T)


common_cols=Reduce(intersect, list(colnames(oligo), colnames(opc), colnames(macro),  colnames(astro_ret), colnames(astro_on)  , colnames(micro), colnames(mg)))

data=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_with_plot/DLG.full.output_ONONH_glia.rds")
mcp=c("MCP1","MCP2","MCP3","MCP4","MCP5")
#mcp=c("MCP5")
for(m in mcp){
oligo_up=oligo[data$MCPs[[m]]$oligo.up,common_cols]
oligo_down=oligo[data$MCPs[[m]]$oligo.down,common_cols]

opc_up=oligo[data$MCPs[[m]]$opc.up,common_cols]
opc_down=oligo[data$MCPs[[m]]$opc.down,common_cols]

macro_up=oligo[data$MCPs[[m]]$macro.up,common_cols]
macro_down=oligo[data$MCPs[[m]]$macro.down,common_cols]

micro_up=oligo[data$MCPs[[m]]$micro.up,common_cols]
micro_down=oligo[data$MCPs[[m]]$micro.down,common_cols]

mg_up=oligo[data$MCPs[[m]]$mg.up,common_cols]
mg_down=oligo[data$MCPs[[m]]$mg.down,common_cols]

astroon_up=oligo[data$MCPs[[m]]$astroon.up,common_cols]
astroon_down=oligo[data$MCPs[[m]]$astroon.down,common_cols]

astroret_up=oligo[data$MCPs[[m]]$astroret.up,common_cols]
astroret_down=oligo[data$MCPs[[m]]$astroret.down,common_cols]

data2=rbind(oligo_up,opc_up, macro_up, micro_up, mg_up, astroon_up, astroret_up, oligo_down, opc_down, macro_down,  micro_down,  mg_down,  astroon_down,  astroret_down)
data1=t(data2)
#DT1=as.matrix(t(scale(t(data1),center=T)))
DT1=as.matrix(scale(data1,center=T))
my_sample_col <- data.frame(celltype =c(rep("OLIGO",dim(oligo_up)[1]), rep("OPC",dim(opc_up)[1]), rep("Macro",dim(macro_up)[1]), rep("Micro",dim(micro_up)[1]),  rep("MG",dim(mg_up)[1]), rep("Astro_ON",dim(astroon_up)[1]), rep("Astro_ret",dim(astroret_up)[1]), rep("OLIGO",dim(oligo_down)[1]),  rep("OPC",dim(opc_down)[1]), rep("Macro",dim(macro_down)[1]), rep("Micro",dim(micro_down)[1]), rep("MG",dim(mg_down)[1]),rep("Astro_ON",dim(astroon_down)[1]), rep("Astro_ret",dim(astroret_down)[1]) ), dir=c(rep("Up",dim(oligo_up)[1]), rep("Up",dim(opc_up)[1]),  rep("Up",dim(macro_up)[1]),   rep("Up",dim(micro_up)[1]), rep("Up",dim(mg_up)[1]),rep("Up",dim(astroon_up)[1]), rep("Up",dim(astroret_up)[1]),  rep("Down",dim(oligo_down)[1]), rep("Down",dim(opc_down)[1]),  rep("Down",dim(macro_down)[1]),rep("Down",dim(micro_down)[1]), rep("Down",dim(mg_down)[1]),rep("Down",dim(astroon_down)[1]), rep("Down",dim(astroret_down)[1]) ))
row.names(my_sample_col) <- colnames(DT1)
ha=HeatmapAnnotation(Celltype=my_sample_col[colnames(DT1),]$celltype,Exp_Dir=my_sample_col[colnames(DT1),]$dir, col=list(Celltype=c("OLIGO"="steelblue","OPC"="pink", "Macro"="darkgreen", "Micro"="yellow", "MG"="lightgreen","Astro_ON"="thistle", "Astro_ret"="grey"), Exp_Dir=c("Down"="lightblue","Up"="pink")),annotation_name_side = "right",annotation_name_gp= gpar(fontsize = 15))
setwd("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_with_plot/")
pdf(paste0("8heatmap_",m,"_glia.pdf"),width=12,height=6)
ht=Heatmap(DT1,row_title = paste0(m,"_glia"),  clustering_distance_rows = "euclidean", show_row_names = FALSE, name = paste0(m,"_glia"),  top_annotation = ha, cluster_columns = F)
draw(ht)
dev.off()
}


#n=match(c
#oligo1=oligo[,(colnames(oligo) %in% colnames(opc))&&(colnames(oligo) %in% colnames(macro))&& (colnames(oligo) %in% colnames(astro)) && (colnames(oligo) %in% colnames(micro)) && (colnames(oligo) %in% colnames(mg))  ]
