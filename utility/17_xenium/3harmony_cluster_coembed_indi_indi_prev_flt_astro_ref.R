library(Seurat)

rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple.rds"
meta1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple_sort.obs.gz"

ref.combined=readRDS(rds) #readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.rds")
meta=read.table(meta1, header=T,sep=",",row.names=1) #read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.obs.gz",header=T,sep="\t",row.names=1)

ref.combined <- AddMetaData(ref.combined, metadata = meta)

###ref.combined = subset(ref.combined, subset = sampleid == args[5])
ref.combined = subset(ref.combined, subset = sampleid %in% c("BCM_22_0047_ONH_RNA","BCM_22_0047_ON_RNA") )

on <- 2200
cells_on <- sample(Cells(ref.combined[,ref.combined@meta.data$subclass2=="ON"]), size = on)

onh=360
cells_onh <- sample(Cells(ref.combined[,ref.combined@meta.data$subclass2=="ONH"]), size = onh)

ret=100
cells_ret <- sample(Cells(ref.combined[,ref.combined@meta.data$subclass2=="retina"]), size = ret)

cells=c(cells_on, cells_onh, cells_ret)

ref=subset(ref.combined, cells = cells)

saveRDS(ref,"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_BCM_22_0047_ON_ONH_RNA_downsample.rds")
