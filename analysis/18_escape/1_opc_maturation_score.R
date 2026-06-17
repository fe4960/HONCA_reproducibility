library(escape)
library(Seurat)
library(BiocParallel)
library(GSVA)
library(msigdbr)
library(dplyr)
name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat"
scRep_example=readRDS(paste0(name, "_simple.rds"))
#GS.hallmark <- getGeneSets(library = "H", subcategory="GO:BP")
#GS.hallmark <- getGeneSets(library = "H",species="Homo sapiens") #, subcategory="GO:BP") #%>% split(x = .$gene_symbol, f = .$gs_name)
param <- BiocParallel::SnowParam(workers = 1)
GS=list(MOL=c("PLP1","MBP","MOG","GRM3","PTGDS","ST18","CERCAM","TMEM144","MOBP","SH3GL3"), OPC=c("PTPRZ1","PDGFRA","VCAN","SOX6","CA10","DCC","NRXN1","TNR"))
scRep_example <- runEscape(input.data=scRep_example, 
                           gene.sets = GS,
			   method = "UCell",
                           groups = 5000, 
                           min.size = 5,
                           new.assay.name = "escape.UCell",
			   BPPARAM = param)

scRep_example <- performNormalization(scRep_example,
                                      assay = "escape.UCell",
                                      gene.sets = GS)

scRep_example <- performNormalization(scRep_example,
                                      assay = "escape.UCell",
                                      gene.sets = GS,
                                      scale.factor = scRep_example$nFeature_RNA)

saveRDS(scRep_example,file=paste0(name, "_mature_escape.rds"))

DefaultAssay(scRep_example)="escape.UCell_normalized"
library(ggplot2)
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
pdf(paste0(name, "_escape_featureplot_UCell_MOL.pdf"))
FeaturePlot(scRep_example, "MOL") +
  scale_color_gradientn(colors = colorblind_vector) +
  theme(plot.title = element_blank())
dev.off()

pdf(paste0(name, "_escape_featureplot_UCell_OPC.pdf"))
FeaturePlot(scRep_example, "OPC") +
  scale_color_gradientn(colors = colorblind_vector) +
  theme(plot.title = element_blank())
dev.off()

