#cell=c("Oligodendrocyte","Astrocyte","Fibroblast","Rod","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","Melanocyte","HC","Cone","RGC","Schwann_cell")
#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas"

args=commandArgs(trailingOnly=T)
cell=read.table(args[1],header=F)$V1
#dir="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas/"
dir=args[2]

prop=function(x){
n=length(x[x<0.75])
s=length(x)
p=n/s
}
for(c in cell){
file=paste0(dir,"/exp_",c,"_cor_Norm_all")

mt=read.table(file,header=T)
mt$a=apply(mt,1,prop)
id=rownames(mt[mt$a>0.65,])

write.table(id, file=paste0(dir,"/exp_",c,"_cor_Norm_all_flt"), quote=F, row.names=F)


}
