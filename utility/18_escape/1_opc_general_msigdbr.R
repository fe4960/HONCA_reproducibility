library(escape)
library(Seurat)
library(BiocParallel)
library(GSVA)
library(msigdbr)
library(dplyr)
#name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat"
args <- commandArgs(trailingOnly = TRUE)
name=args[1]

scRep_example=readRDS(paste0(name, ".rds"))

#scRep_example=readRDS(paste0(name, "_simple.rds"))
#GS.hallmark <- getGeneSets(library = "H", subcategory="GO:BP")
GS.hallmark <- msigdbr(species  = "Homo sapiens",  
                       category = "H")

#GS.hallmark <- getGeneSets(library = "H",species="Homo sapiens") #, subcategory="GO:BP") #%>% split(x = .$gene_symbol, f = .$gs_name)
#param <- BiocParallel::SnowParam(workers = 2)
param <- BiocParallel::SnowParam(workers = 1)

scRep_example <- runEscape(input.data=scRep_example, 
                           gene.sets = GS.hallmark,
			   method = "ssGSEA",
                           groups = 5000, 
                           min.size = 5,
                           new.assay.name = "escape.ssGSEA",
#			   msigdbGeneSets = c("H"),
			   BPPARAM = param)

scRep_example <- performNormalization(scRep_example,
		#		      method= "ssGSEA",
                                      assay = "escape.ssGSEA",
		#		      min.size=5,
                                      gene.sets = GS.hallmark)
#				      msigdbGeneSets = c("H"),
#				      )

scRep_example <- performNormalization(scRep_example,
                                      assay = "escape.ssGSEA",
		#		      method= "ssGSEA",
		#		      min.size=5,
                                      gene.sets = GS.hallmark,
#				      msigdbGeneSets = c("H"),
                                      scale.factor = scRep_example$nFeature_RNA)

saveRDS(scRep_example,file=paste0(name, "_escape_msigdbr.rds"))

library(ggplot2)
library(RColorBrewer)
pdf(paste0(name, "_escape_heatmap_msigdbr.pdf"))
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = "all", #rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.ssGSEA",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 
dev.off()


pdf(paste0(name, "_escape_heatmap_norm_msigdbr.pdf"))
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = "all", #rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.ssGSEA_normalized",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 
dev.off()





all.markers <- FindAllMarkers(scRep_example,
			      group.by="subclass",
                              assay = "escape.ssGSEA_normalized",
                              min.pct = 0,
                              logfc.threshold = 0)

df=data.frame(all.markers)

write.table(df,file=paste0(name, "_escape_DEG_msigdbr.txt"),sep="\t",quote=F)

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
DefaultAssay(scRep_example)="escape.ssGSEA_normalized"


#pdf(paste0(name, "_escape_featureplot_G2M.pdf"))
#FeaturePlot(scRep_example, "HALLMARK-G2M-CHECKPOINT") +
#  scale_color_gradientn(colors = colorblind_vector) +
#  theme(plot.title = element_blank())
#dev.off()

pdf(paste0(name, "_escape_year_msigdbr.pdf"))
RidgePlot(
  scRep_example,
  "age_year",
  group.by = "subclass",
)
dev.off()
