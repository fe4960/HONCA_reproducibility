library(Seurat)
library(DIALOGUE)
obj=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/ONONH_majorclass_40000.rds")
#obj1=readRDS("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/ONONH_Astrocyte_celltype_20000.rds")

DIALOGUE_make.cell.type.seurat<- function(obj, cell.subtypes, name, type){
  # Given a Seurat object with "cell.subtypes" and "samples" information provided in the meta.data field
  
#if(is.null(obj@meta.data$cell.subtypes)){return("Cell subtypes information is missing")}
#if(is.null(obj@meta.data$samples)){return("Sample information is missing")}
  
data=subset(x = obj,cells = rownames(obj@meta.data)[is.element(obj@meta.data[[type]],cell.subtypes)])

data<- NormalizeData(data,normalization.method = "LogNormalize", scale.factor = 10000)
data<-FindVariableFeatures(data,selection.method = "vst", nfeatures = 2000)
data<-ScaleData(data)
data<-RunPCA(data,npcs = 30) 
  
#tpm<- as.matrix(sub_obj@assays$RNA@data

tpm=as.matrix(data@assays$RNA@data)
X=data@reductions$pca@cell.embeddings
assertthat::are_equal(rownames(X),colnames(tpm))
n_cells=ncol(tpm)
samples=data@meta.data$donor
n_count<- as.vector(scale(log(data@meta.data$nCount_RNA)))
#gender=as.character(data@meta.data$gender)
#race=as.character(data@meta.data$race)
#tissue=as.character(data@meta.data$tissue)
ageyear=as.numeric(data@meta.data$age_year)
df=data.frame(age=ageyear,agerange="60younger")
df[df$age>=60,]$agerange="60orMore"

metadata<- data.frame(nFeatures = data@meta.data$nFeature_RNA,
                        samples = samples,
                        cell.subtypes = as.character(data@meta.data[[type]]),
#			cell.celltypes= as.character(data@meta.data$celltype),
			gender=as.character(data@meta.data$gender),
			race=as.character(data@meta.data$race),
			tissue=as.character(data@meta.data$tissue),
			ageyear=as.numeric(data@meta.data$age_year),
			age=as.character(data@meta.data$age),
			donor=as.character(data@meta.data$donor),
			agerange=as.character(df$agerange))
r<-make.cell.type(name=name ,tpm,samples,X,metadata, cellQ=n_count )
return(r)
}

cell=c("Astrocyte","Endothelial_cell","Fibroblast","MG","Macrophage","Melanocyte","Microglia","Mural_cell","Oligodendrocyte","Oligodendrocyte_precursor_cell","RPE","RGC")

#rA=c()

#for(cell.subtypes in cell){
#cell.subtypes=c("Astrocyte","MG","Macrophage","Microglia")
#r_astro=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Astrocyte"), name="astro_retina")

r_micro=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Microglia"), name="micro", type="majorclass")
#r_macro=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Macrophage"), name="macro", type="majorclass")
r_oligo=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Oligodendrocyte"), name="oligo", type="majorclass")
r_opc=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Oligodendrocyte_precursor_cell"), name="opc",type="majorclass")
#r_mg=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("MG"), name="MG",type="majorclass")
r_astro=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Astrocyte"), name="astro",type="majorclass")

#r_atro_ret_nlgn=DIALOGUE_make.cell.type.seurat(obj=obj1, cell.subtypes=c("Astro_retina_NLGN1+"), name="astro_ret_nlgn1",type="celltype")
#r_astro_ret=DIALOGUE_make.cell.type.seurat(obj=obj1, cell.subtypes=c("Astro_ONH","Astro_retina","Astro_retina_NLGN1+"), name="astro_ret",type="celltype")
#r_astro_onh=DIALOGUE_make.cell.type.seurat(obj=obj1, cell.subtypes=c("Astro_ONH"), name="astro_onh",type="celltype")
#r_astro_ononh=DIALOGUE_make.cell.type.seurat(obj=obj1, cell.subtypes=c("Astro_ONONH"), name="astro_ononh",type="celltype")
#r_astro_on=DIALOGUE_make.cell.type.seurat(obj=obj1, cell.subtypes=c("Astro_ONONH","Astro_ON","Astro_ON_AFF3+"), name="astro_on",type="celltype")
#r_atro_on_aff3=DIALOGUE_make.cell.type.seurat(obj=obj1, cell.subtypes=c("Astro_ON_AFF3+"), name="astro_on_aff3",type="celltype")

#r_rgc=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("RGC"), name="RGC")

#r_oligo=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Oligodendrocyte","Oligodendrocyte_precursor_cell"), name="oligo")
#r_glia_all=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Astrocyte","MG","Macrophage","Microglia","Oligodendrocyte","Oligodendrocyte_precursor_cell"), name="glia_all")
#r_glia_rgc=DIALOGUE_make.cell.type.seurat(obj=obj, cell.subtypes=c("Astrocyte","MG","Macrophage","Microglia","Oligodendrocyte","Oligodendrocyte_precursor_cell","RGC"), name="glia_rgc")

#rA<- list(r_astro_on=r_astro_on, r_astro_ononh = r_astro_ononh, r_astro_onh = r_astro_onh  , r_astro_ret = r_astro_ret, r_micro = r_micro, r_macro=r_macro, r_oligo=r_oligo,r_opc=r_opc, r_mg=r_mg)
rA<- list(r_astro=r_astro, r_micro = r_micro,  r_oligo=r_oligo,r_opc=r_opc) #, r_mg=r_mg)

#rA=c(r_tmp,rA)
#}

######param<-DLG.get.param(k = 2, seed1 = 7, results.dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_2k_Age/",  n.genes = 100, conf = "cellQ", spatial.flag = F,plot.flag = F , pheno = "agerange")


####param<-DLG.get.param(k = 3, seed1 = 7, results.dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_3k_Age/",  n.genes = 100, conf = "cellQ", spatial.flag = F,plot.flag = T ,  pheno = "ageyear")
param<-DLG.get.param(k = 3, seed1 = 7, results.dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_3k_AgeRange/",  n.genes = 100, conf = "cellQ", spatial.flag = F,plot.flag = T ,  pheno = "agerange")


R <- DIALOGUE.run(rA = rA, # list of cell.type objects
                    main = "ONONH_glia_only",
                    param = param)
saveRDS(R,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_3k_AgeRange/glia_only_on_noPlot.rds")

#saveRDS(R,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_3k_Age/glia_only_on_noPlot.rds")

#saveRDS(R,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/glia_only_on_noPlot_2k_Age/glia_only_on_noPlot.rds")
#pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/DIALOGUE/plot.pdf",width=15,height=20)
#  par(mfrow=c(3,2),oma = c(8, 1, 0, 5),xpd = T)
  # MCP1 marks CD8 T cells and DCs
#  boxplot(-R$scores$r_glia$MCP1~R$scores$r_glia$cell.subtypes,xlab = "",ylab = "MCP1",las=2)
#  boxplot(-R$scores$r_stroma$MCP1~R$scores$r_stroma$cell.subtypes,xlab = "",ylab = "MCP1",las=2)
#  boxplot(-R$scores$r_oligo$MCP1~R$scores$r_oligo$cell.subtypes,xlab = "",ylab = "MCP1",las=2)
#  boxplot(-R$scores$r_glia_all$MCP1~R$scores$r_glia_all$cell.subtypes,xlab = "",ylab = "MCP1",las=2)
#  boxplot(-R$scores$r_glia_rgc$MCP1~R$scores$r_glia_rgc$cell.subtypes,xlab = "",ylab = "MCP1",las=2)

#  return(R)
#dev.off()

