library(ggplot2)
library(patchwork)

#files=c("Oligodendrocyte_precursor_cell_subclass_seurat_simple","Oligodendrocyte_subclass_seurat_simple" ) #RPEchoroid_scRNA_majorclass_5k
#files=c("Fibro_TM","Fibro_CB")
files=c("Fibro_TM", "cb_fibro_subtype_02_rawcount_anno_clean_simple" )
value1=NULL
for(file in files){
value=read.table(paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list1_summary_list_new_MAGMA_",file),header=T)

value$Trait=gsub("POAG_GERA_UKB2018","POAG_GERA",value$Trait)
value$Trait=gsub("AMD_NG2016","AMD",value$Trait)
value$Trait=gsub("Khor_PACG2016","PACG",value$Trait)
value$Trait=gsub("Myopia2018NG","Myopia",value$Trait)
value$Trait=gsub("CA_DA_caucasians_2017","POAG_CA",value$Trait)
value$Trait=gsub("DA_caucasians_2017","POAG_DA",value$Trait)
value$Trait=gsub("IOP_caucasians_2017","POAG_IOP",value$Trait)
value$Trait=gsub("VCDR_caucasians_2017","POAG_VCDR",value$Trait)
value$Trait=gsub("UKBBc_IOP_2018Khawaja","IOP",value$Trait)
value$Trait=gsub("POAG_Gharahkhani2021","POAG",value$Trait)
value$Trait=gsub("Refracive_Error_Hysi2020","Refractive_Error",value$Trait)
value$Trait=gsub("Bone_mineral_density","Bone_mineral_density",value$Trait)
value$Trait=gsub("WBC","WBC",value$Trait)
value$Trait=gsub("UKBB_BMI","BMI",value$Trait)
value$Cell_type=gsub("RBC","BCR",value$Cell_type)
value$Cell_type=gsub("Astrocyte","Astro",value$Cell_type)
value$Cell_type=gsub("Microglia","Micro",value$Cell_type)
value$Cell_type=gsub("Microglia","Micro",value$Cell_type)
value$Cell_type=gsub("Oligodendrocyte_precursor_cell","OPC",value$Cell_type)

value$Trait=gsub("UKBB_IOP","IOP",value$Trait)

value=value[value$Trait%in%c("PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","IOP","Bone_mineral_density"),]
#value=value[value$Trait%in%c("Disorders_of_choroid_and_retina","Retinal_detachments_and_breaks","Retinal_problem","AMD","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","Diabetic_retinopathy","ONL_thickness","IS_thickness","OS_thickness","Bone_mineral_density", "UKBB_IOP"),]

###value$p.adj=p.adjust(value$FDR,method="BH")

value$p.adj=p.adjust(value$P.val,method="BH")
#value$color="darkseagreen1"
value$color="grey"
value[value$FDR<0.1,]$color="mediumpurple"

#value[value$p.adj<0.05,]$color="darkgreen"
#value[value$P.val<0.05,]$color="mediumpurple"
#value[value$P.val<0.05,]$color="#8DA0CB"
#value[value$FDR<0.05,]$color="purple"

value$ct=file
value1=rbind(value,value1)
}



value=value1
p <- ggplot(data=value, aes(factor(Cell_type),factor(Trait,levels=c("AMD","PACG","POAG","POAG_CA","POAG_DA","POAG_VCDR","POAG_IOP","IOP","UKBB_IOP","Diabetic_retinopathy","ONL_thickness","IS_thickness","OS_thickness","Disorders_of_choroid_and_retina","Retinal_detachments_and_breaks","Retinal_problem","Bone_mineral_density"))))
##p_dar=p + geom_point(aes(size=-log(value$FDR,base=10)),color=value$color)+theme(

p_dar=p + geom_point(aes(size=-log(value$P.val,base=10)),color=value$color)+theme(
    axis.line=element_line(colour="black")
    ,panel.background = element_rect(fill="transparent") # bg of the panel
    , plot.background = element_blank() # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    #, legend.position = "none"
    , text=element_text(size=28),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
####  )+xlab("")+ylab("")+  guides(color = guide_colorbar(order = 1),fill = guide_legend(order = 1))+scale_size_continuous(name=expression(paste('-',log[10],'FDR')),range=c(1,15))+facet_wrap(~ct, scales = "free_x")

   # axis.text.x = element_text(size=10)  # get rid of legend bg
  )+xlab("")+ylab("")+  guides(color = guide_colorbar(order = 1),fill = guide_legend(order = 1))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(1,15))+facet_wrap(~ct, scales = "free_x")

####  guides(color = guide_colorbar(order = 1),fill = guide_legend(order = 1))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(1,15))+facet_wrap(~ct, scales = "free_x")
###pdf(paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list_summary_list_new_MAGMA_oligo_opc_comb_FDR.pdf"),width=11)


pdf(paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list_summary_list_new_MAGMA_TMCB_Fibro1.pdf"),width=11)

#pdf(paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list_summary_list_new_MAGMA_oligo_opc_comb.pdf"),width=11)
 print(p_dar)
 dev.off()
#### write.table(value,file=paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list_summary_list_new_MAGMA_oligo_opc_padj_FDR"),quote=F,sep="\t")
write.table(value,file=paste0("/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list_summary_list_new_MAGMA_TMCB_Fibro_padj1"),quote=F,sep="\t")

