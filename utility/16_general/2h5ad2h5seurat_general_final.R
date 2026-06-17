library(Seurat)
library(SeuratDisk)
#library(rhdf5)
args <- commandArgs(trailingOnly = TRUE)
sam=args[1]
dir=args[2]
#sam1=args[3]
#dir="/dfs3b/ruic20_lab/junw42/human_ret_anc/data/snRNA/"
Convert(paste0(dir ,sam,".h5ad"), dest = "h5seurat", overwrite = TRUE)

seRNA0 <- LoadH5Seurat(paste0(dir,sam,".h5seurat"),misc = FALSE,meta.data = FALSE)

meta.data=read.table(paste0(dir, sam ,".obs.gz"),header=T,sep=",",row.names=1)

#meta.data=read.table(paste0(dir, sam1 ,".obs.gz"),header=T,sep=",",row.names=1)

#rownames(meta.data)=meta.data$barcode

seRNA=AddMetaData(object=seRNA0, metadata=meta.data)

saveRDS(seRNA,file=paste0(dir, sam,".rds"))

