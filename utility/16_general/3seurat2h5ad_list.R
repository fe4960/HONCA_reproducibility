
library(Seurat)
library(SeuratDisk)
#library(rhdf5)
args <- commandArgs(trailingOnly = TRUE)
sam=args[1]
dir=args[2]
#rds=args[3]
#sam1=args[3]
#dir="/dfs3b/ruic20_lab/junw42/human_ret_anc/data/snRNA/"
#Convert(paste0(dir ,sam,".h5ad"), dest = "h5seurat", overwrite = TRUE)

obj=readRDS(paste0(dir,"/", sam,".rds"))

for(i in length(obj)){
obj1=obj[[i]]
filename=paste0(dir,"/", sam,"_",i,".h5Seurat") #output_PATH, target_tissue, "_coembed.combined_inte.h5Seurat")
SaveH5Seurat(obj1, filename = filename,overwrite = TRUE)
Convert(filename, dest = "h5ad",overwrite = TRUE)
}
