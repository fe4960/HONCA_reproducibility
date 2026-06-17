library(ggplot2)


generate_txt=function(outdir, bname){
file=paste0(outdir, "/", bname, "_race.csv")
data=read.table(file, header=T, sep="\t")
p=ggplot(data=data,aes(x=subclass,y=nCount_RNA,fill=race))+geom_bar(stat = "identity",position="fill")+labs(x = "Ancestry", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_brewer(palette = "Paired") #values = c("yellow", "purple"))
pdf(paste0(outdir, bname,"_race_df.pdf"))
print(p)
dev.off()

file=paste0(outdir, "/", bname, "_gender.csv")
data=read.table(file, header=T, sep="\t")
p=ggplot(data=data,aes(x=subclass,y=nCount_RNA,fill=gender))+geom_bar(stat = "identity",position="fill")+labs(x = "Sex", y = "Cell proportion") + theme_minimal()+theme(axis.text.x = element_text(angle = 65, hjust = 1), text=element_text(size = 25), plot.margin = margin(t = 10, r = 10, b = 25, l = 10))+scale_fill_manual(values = c("yellow", "purple"))
pdf(paste0(outdir, bname,"_gender_df.pdf"))
print(p)
dev.off()
}


outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
bname="Astrocyte_subclass_new5_clean"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
bname="Oligodendrocyte_precursor_cell_subclass_seurat_cycling"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
bname="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location"
generate_txt(outdir, bname)

outdir='/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/'
bname='Endothelial_cell_subclass_sb_seurat_rmRPE'
generate_txt(outdir, bname)


outdir='/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/'
bname="Mural_cell_subclass_rmRPE"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
bname="Oligodendrocyte_subclass_seurat"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
bname="Macrophage_subclass_sb_clean_rmRPE"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
bname="Microglia_subclass_sb_clean"
generate_txt(outdir, bname)


