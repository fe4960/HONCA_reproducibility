# Usage: Rscript /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_250317/DRG_Xenium_202503/9a_coembed_harmony.R xenium.combined neuron/nonneuron/all

library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)


od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
#kw="PP_" #args[4]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
id="5" #args[6]
#getwd()
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue1, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/")
dir.create(output_PATH, recursive = TRUE)
# read xenium data
#path=paste0("/dfs3b/ruic20_lab/rawdata/Xenium/",path0)

#xenium.combined <- args[1]
#parts <- args[2]

#============================================================================================================
# co-embedding - Harmony


library(harmony)
library(RColorBrewer)




#merged_obj=readRDS(paste0(output_PATH, target_tissue, "_8_coembed.combined_cutoff20_ref_xenium_no_inte_harmony.rds"))
	


xenium.combined=readRDS(paste0(output_PATH, target_tissue, "_8_coembed.combined_inte_harmony_w1", ".rds"))

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_8_xenium_data_list_comb_cutoff20.rds"))

cut=20

#xenium.combined=merged_obj   #readRDS(paste0(output_PATH, target_tissue, "coembed.combined_inte_harmony_w1", ".rds"))   #readRDS( paste0(output_PATH, target_tissue, "_xenium.combined_cutoff",cut,".rds"))

#id=unique(xenium.combined@meta.data$slide_id)

library(ggplot2)
#for(i in 1:(length(id)-1)){
#for(i in 1:(length(id)){
for(i in 1:length(xenium_data_list)){

xenium_obj=xenium_data_list[[i]]

slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

merged_obj1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]


#meta=x1@meta.data #merged_obj@meta.data[merged_obj@meta.data$slide_id==id[i],]

xenium_obj@meta.data$harmony_anno="other"
xenium_obj@meta.data[merged_obj1@meta.data$barcode,]$harmony_anno=merged_obj1@meta.data$harmony_anno

xenium_obj@meta.data[xenium_obj@meta.data$harmony_anno=="Mast_cell",]$harmony_anno="Immune_cell"

DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

#paired_more = c("azure4", brewer.pal(8, "Set2"), "darkviolet", brewer.pal(12, "Set3"), brewer.pal(8, "Set1"))
paired_more = c("azure4", brewer.pal(8, "Set2"), "darkviolet", brewer.pal(12, "Set3"), brewer.pal(8, "Set1"))  # w1_2


#paired_more = c("azure4", brewer.pal(8, "Set1"), brewer.pal(12, "Set3"),  "darkviolet",   brewer.pal(8, "Set2"))  # w1_3

#paired_more <- c(
#  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
#  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
#  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
#  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
#  "#969696", "#DD1C77", "#6BAED6", "#74C476", "#9E9AC8",
#  "#FD8D3C", "#BDBDBD", "#FDD0A2", "#AEC7E8"
#) # w1_4


#paired_more = c("azure4", brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(8, "Set1"))

celltypes <- sort(unique(xenium_obj$harmony_anno))

#########cols_named <- setNames(paired_more[1:length(celltypes)], celltypes)


#cols_named=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_8_all_annot_majorclass_ONH01_w1_2_color1", header=T, sep="\t", row.names=1)
# Save named colors to a tab-delimited text file
# light blue "#80B1D3",
cols_named <- c(
"AC" =     "azure4",
"Astrocyte"   =   "#A6D854", #"#66C2A5",
"BC"   =   "#FC8D62",
"Cone" =    "#8DA0CB",
"Endothelial_cell" =        "#E78AC3",
"Fibroblast" =  "#00E5FF",  # "#A6D854",
"HC" =      "#FFD92F",
"Immune_cell" =     "#E5C494",
"Macrophage"  =    "#B3B3B3",
"Mast_cell"  = "#3399FF", # "#A6D854",  
"Melanocyte"   =   "#8DD3C7",
"MG"="#FFFFB3",
"Microglia"="#BEBADA",
"Mural_cell" = "#FDE725", # "#FB8072",
"Oligodendrocyte"="darkviolet",
"Oligodendrocyte_precursor_cell" = "#FDB462",
"RGC" = "#B3DE69",
"Rod" =     "#FCCDE5",
"RPE" = "#D9D9D9",
"Schwann_cell" = "#B36A1A" #"#BC80BD"
)

##########
#########
#cols_named= c( "Melanocyte"   =   "#8DD3C7",
#  "Oligodendrocyte_precursor_cell" = "#FDB462",
#  "Mast_cell"  =  "#A6D854",
#  "Schwann_cell" ="azure4",
#  "AC" = "#E64B65",
#  "BC" = "#00E5FF",
#  "Cone" = "#00A087",
#  "HC" = "#91D1C2",
#  "Rod"  =  "#00A8C6",
#  "MG" =   "#4FC3F7",
#  "RPE" = "#B39DDB",
#  "Astrocyte" = "#7AD151",
#  "Oligodendrocyte" = "#4CAF50",
#  "Fibroblast" = "#FDE725",
#  "Endothelial_cell" = "#B8DE29",
#  "Mural_cell" = "#F39B7F",
#  "Microglia" = "#D55E00",
#  "Macrophage" = "#B36A1A",
#  "Endothelial_cell" = "#CC79A7",
#  "Immune_cell" = "#3399FF",
#  "RPE" = "#FF69B4")









#  "B_Cell" = "#E64B65",
#  "NK/T cell" = "#00E5FF",
#  "Macrophage" = "#00A087",
#  "Mast_Cell" = "#91D1C2",
#  "Schwann-my"  =  "#00A8C6",
#  "Schwann-nmy" =   "#4FC3F7",
#  "Ciliary_Muscle" = "#B39DDB",
#  "CB_PCE_1" = "#7AD151",
#  "CB_PCE_2" = "#4CAF50",
#  "CB_NPCE_1" = "#FDE725",
#  "CB_NPCE_2" = "#B8DE29",
#  "CB_Fibro_BMP5+" = "#F39B7F",
#  "CB_Fibro_CNTN4+" = "#D55E00",
#  "CB_Fibro_SMOC2+" = "#B36A1A",
#  "Uveal_Melanocyte" = "#CC79A7",
#  "Vascular_Endothelium" = "#3399FF",
#  "Pericyte1" = "#FF69B4")

write.table(
  data.frame(
    celltype = names(cols_named),
    color = cols_named
  ),
  file = paste0(output_PATH, target_tissue, "_8_all_annot_majorclass_", slide_id, "_w1_7_color"),

	    #  file = paste0(output_PATH, target_tissue, "_8_all_annot_majorclass_", slide_id, "_w1_5_color"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#pdf(paste0(output_PATH, target_tissue, "_8_all_annot_majorclass_", slide_id, "_w1_3.pdf"), width = 50, height = 35)
#pdf(paste0(output_PATH, target_tissue, "_8_all_annot_majorclass_", slide_id, "_w1_4.pdf"), width = 50, height = 35)

pdf(paste0(output_PATH, target_tissue, "_8_all_annot_majorclass_", slide_id, "_w1_7.pdf"), width = 50, height = 35)

#pdf(paste0(output_PATH, target_tissue, "_8_all_annot_majorclass_", slide_id, "_w1_5.pdf"), width = 50, height = 35)
p=ImageDimPlot(object = xenium_obj, group.by = "harmony_anno", border.size = 0,   cols = cols_named, border.color=NA, axes=TRUE, dark.background=TRUE )
#"polychrome"
#+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
print(p)
dev.off()

}

#####

#========

