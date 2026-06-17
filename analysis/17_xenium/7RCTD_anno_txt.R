library(Seurat)
library(spacexr)


od="20250620_PP_comb" #args[1]
target_tissue="ONONH_comb" #args[2]
target_tissue1="ONONH"

path0="20250620__203138__Human_70y_ASRetina_20250620" #args[3]
cb_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/" #args[5]
wd=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/",od)
setwd(wd)
output_PATH1 <- paste0(wd, "/",target_tissue, "/")
output_PATH <- paste0(wd, "/",target_tissue, "/spacexr/")
dir.create(output_PATH, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

label=args[1]
clu=args[2]
idx=args[3]
sl=c("PP_1","PP_2","ONH_1","ONH_2","ONH_3")
for(t in sl){
fn=paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron",idx,".rds")
fn1=paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron",idx,".txt")
if(idx == "n"){
fn=paste0(output_PATH, label,"_", t, "_RCTD_results_nonneuron.rds")

#fn=paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron.rds")
fn1=paste0(output_PATH, label,"_", t, "_",clu,"_RCTD_results_nonneuron.txt")

}
RCTD=readRDS(fn)
data=RCTD@results$results_df

write.table(data, file=fn1, sep="\t", quote=F)
}
