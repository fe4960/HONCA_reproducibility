library(Seurat)
library(ggplot2)
library(dplyr)
library(spacexr)
#library(arrow)

options(future.globals.maxSize = 50 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)
od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
############
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH <- paste0(wd, "/",target_tissue, "/")

cut=as.numeric(args[2])

########
#####process xenium data

label=args[1]

fc=as.numeric(args[3])
clu=args[4]
sc=args[5]

xenium_data_list=readRDS(paste0(output_PATH, target_tissue,"_xenium_data_list_comb_cutoff20.rds"))


xenium.combined=readRDS(paste0(output_PATH, target_tissue,"_",label, "_xenium.combined_cutoff",cut,"_fc_",fc,".rds"))


id=unique(xenium.combined@meta.data$slide_id)

x1_meta_all=NULL

for(i in 1:length(id)){


xenium_obj=xenium_data_list[[i]]
slide_id=unique(xenium_data_list[[i]]@meta.data$slide_id)

x1=xenium.combined[,xenium.combined@meta.data$slide_id==slide_id]

xenium_obj=xenium_obj[,x1@meta.data$barcode]


#xenium_obj@meta.data[,clu]=x1@meta.data[,clu]


RCTD=readRDS(paste0(output_PATH, "/spacexr/",label,"_", slide_id,"_",sc, "_RCTD_results_nonneuron.rds"))

#xenium_obj@meta.data[,clu]=x1@meta.data[,clu]
x1@meta.data[,"rctd1"]=RCTD@results$results_df[x1@meta.data$barcode,]$first_type
x1=x1[,x1@meta.data[,"rctd1"]!="NA"]
#DefaultBoundary(xenium_obj[["fov"]]) <-  "segmentation"

#pdf(paste0(output_PATH, target_tissue, "_",label,"_annot_subclass_",id[i],"_cutoff",cut,"_fc_",fc,"_",clu,"_only.pdf"), width = 50, height = 35)
#p=ImageDimPlot(object = xenium_obj, group.by = clu, border.size = 0,  axes=TRUE, dark.background=TRUE, cols = "polychrome")
#print(p)
#dev.off()



x1@meta.data$y_coor=x1@images$fov$centroids$y


if((i==2) || (i==5)){
x1@meta.data$y_coor=-x1@images$fov$centroids$y
}

#library(ggplot2)
############



x1@meta.data$bin_y <- cut(
  x1@meta.data$y_coor,
  breaks = 50,         # number of bins
  labels = FALSE,      # return integers 1..10
  include.lowest = TRUE
)

x1_meta=x1@meta.data

x1_meta_all=rbind(x1_meta_all,x1_meta)

}



library(rlang)

library(dplyr)



df <- x1_meta_all %>%
  group_by(bin_y, rctd1   ) %>%

#  group_by(bin_y, !!sym(clu)  ) %>%
  summarise(
    n = n()
  ) %>%
  group_by(bin_y) %>%
  mutate(prop = n / sum(n))


fn=paste0(output_PATH, target_tissue,"_",label,"_cutoff",cut,"_fc_",fc,"_",clu,"_cell_cluster_dist_y_comb_rmNA.pdf")

#fn=paste0(output_PATH, target_tissue,"_",label,"_cutoff",cut,"_fc_",fc,"_",clu,"_cell_cluster_dist_y_comb.pdf")
#p=ggplot(df,aes(x=bin_y,y=prop,fill= !!sym(clu)  ))+geom_bar(stat = "identity") +
p=ggplot(df,aes(x=bin_y, y=prop, fill=rctd1  ))+geom_bar(stat = "identity") +
ylab("Proportion") + theme_minimal() + scale_fill_brewer(palette="Set3")
pdf(fn)
print(p)
dev.off()

#########


#for(i in 1:length(id)){

#}




#}
#summary(lm(df[df$xenium_clusters_p2==0,]$prop~df[df$xenium_clusters_p2==0,]$bin))
