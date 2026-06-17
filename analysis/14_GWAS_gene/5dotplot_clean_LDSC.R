library(qvalue)
library(ggplot2)


value_all=NULL

files=c("major", "peak")
for(file in files){
#value=read.table(paste0("~/Documents/human_meta/file_list_new_control_",file,"_summary_list"),header=T)
value=read.table(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/14_GWAS_gene//file_list_new_control_",file,"_summary_list"), header=T, sep="\t")
#value$Trait=gsub("Advanced_age-related_macular_degeneration","Age-related_macular_degeneration",value$Trait)
#value$Trait=gsub("Heel_bone_mineral_density","Bone_mineral_density",value$Trait)

#value$qval=qvalue(value$Pval)$qvalues
#value=value[value$Cell_type!="Endo"&value$Cell_type!="Mic",]
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
value$Trait=gsub("Refracive_Error_Hysi2020","Refracive_Error",value$Trait)
value$Trait=gsub("Bone_mineral_density","Bone_mineral_density",value$Trait)
value$Trait=gsub("WBC","WBC",value$Trait)
value$Trait=gsub("UKBB_BMI","BMI",value$Trait)
#value$Cell_type=gsub("RBC","BCR",value$Cell_type)
#value$Cell_type=gsub("Astrocyte","Astro",value$Cell_type)
#value$Cell_type=gsub("Microglia","Micro",value$Cell_type)
#value=value[value$Cell_type!="Endo"&value$Cell_type!="Micro"&value$Cell_type!="RPE",]
#value=value[value$Trait%in%c("Disorders_of_choroid_and_retina","Retinal_detachments_and_breaks","Retinal_problem","AMD","Refracive_Error","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","Bone_mineral_density","Diabetic_retinopathy","ONL_thickness","IS_thickness","OS_thickness"),]
value=value[value$Trait%in%c("AMD","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","Bone_mineral_density","ONL_thickness","IS_thickness","OS_thickness"),]
value=value[value$Cell_type!="MG",]
tr=unique(value$Trait)
for(t in tr){
print(t)
tr1=value[value$Trait==t,]
#tr1$qval=qvalue(as.numeric(tr1$P.val))$qvalue
tr1$FDR=p.adjust(tr1$P.val,method="BH")
tr1$color="grey"
tr1[tr1$FDR<0.1,"color"]="gold"
tr1$type=file
value_all=rbind(value_all,tr1)
}
}
colnames(value_all) = colnames(tr1)
#value$p.adj=p.adjust(value$P.val,method="BH")
#value$qval=qvalue(value$P.val)$qvalue
#value$color="grey"
#value[value$qval<0.1,]$color="orange"
#p <- ggplot(value, aes(factor(Cell_type, levels=c("Rod","Cone","BC","HC","AC","RGC","MG","Astro")),factor(Trait,levels=c("Disorders_of_choroid_and_retina","Retinal_detachments_and_breaks","Retinal_problem","AMD","Refracive_Error","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","Diabetic_retinopathy","ONL_thickness","IS_thickness","OS_thickness","Bone_mineral_density"))))
p <- ggplot(value_all, aes(x=Cell_type,y=factor(Trait,levels=c("AMD","PACG","POAG_CA","POAG_DA","POAG_VCDR","IOP","POAG","ONL_thickness","IS_thickness","OS_thickness","Bone_mineral_density"))))



p_dar=p + geom_point(aes(size=-log(value_all$P.val,base=10)),color=value_all$color)+theme(
    axis.line=element_line(colour="black")
    ,panel.background = element_rect(fill="transparent") # bg of the panel
    , plot.background = element_blank() # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    #, legend.position = "none"
    , text=element_text(size=25),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
   # axis.text.x = element_text(size=10)  # get rid of legend bg
  )+xlab("")+ylab("")+ggtitle(paste0("snATAC-seq ",file))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(1,15))+facet_grid(~type, scales="free_x", space="free_x") #+scale_color_gradient(name=expression(paste('-',log[10],'P')),low = "grey", high = "blue")+scale_color_continuous(name=expression(paste(Coef,'*',10^6)))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(1,15))
file="all"
 pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/14_GWAS_gene/file_list_new_control_",file,"_new_summary_list2.pdf"),width=12, height=8)

# pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/14_GWAS_gene/file_list_new_control_",file,"_new_summary_list.pdf"),width=9, height=8)
 print(p_dar)
 dev.off()
 write.table(value_all,file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/14_GWAS_gene/file_list_new_control_",file,"_summary_list_padj2"),quote=F,sep="\t")

 # write.table(value_all,file=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/14_GWAS_gene/file_list_new_control_",file,"_summary_list_padj"),quote=F,sep="\t")
# }


