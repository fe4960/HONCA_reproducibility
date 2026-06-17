library(escape)
library(Seurat)
library(BiocParallel)
library(GSVA)
library(msigdbr)
library(dplyr)
name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat"
scRep_example=readRDS(paste0(name, "_simple.rds"))
#GS.hallmark <- getGeneSets(library = "H", subcategory="GO:BP")
GS.hallmark <- getGeneSets(library = "H",species="Homo sapiens") #, subcategory="GO:BP") #%>% split(x = .$gene_symbol, f = .$gs_name)
param <- BiocParallel::SnowParam(workers = 2)
#msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# 3. Split into list format (what GSVA and escape expect)
#GS.hallmark <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

#enrichment.scores <- escape.matrix(scRep_example,
#                                   gene.sets = GS.hallmark,
#                                   groups = 5000,
#                                   min.size = 5,
#                                   BPPARAM = SnowParam(workers = 4))


#scRep_example <- performNormalization(scRep_example,
#                                      assay = "escape.ssGSEA",
#                                      gene.sets = GS.hallmark)

#BPPARAM = SnowParam(workers = 4)

#scRep_example <- runEscape(scRep_example,
#                           method = "ssGSEA",
#                           gene.sets = GS.hallmark,
#                           groups = 5000,
#                           min.size = 5,
#                           new.assay.name = "escape.ssGSEA")

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

saveRDS(scRep_example,file=paste0(name, "_escape.rds"))

library(ggplot2)
library(RColorBrewer)
pdf(paste0(name, "_escape_heatmap.pdf"))
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.ssGSEA",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 
dev.off()


pdf(paste0(name, "_escape_heatmap_norm.pdf"))
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.ssGSEA_normalized",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 
dev.off()




#palette = "Spectral"

pdf(paste0(name, "_escape_ridge_angio.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-ANGIOGENESIS")
dev.off()


pdf(paste0(name, "_escape_ridge_interferon.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-INTERFERON-ALPHA-RESPONSE")
dev.off()

pdf(paste0(name, "_escape_ridge_complement.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-COMPLEMENT")
dev.off()

pdf(paste0(name, "_escape_ridge_G2M.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-G2M-CHECKPOINT")
dev.off()

pdf(paste0(name, "_escape_ridge_MYC_V2.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-MYC-TARGETS-V2")
dev.off()

pdf(paste0(name, "_escape_ridge_NOTCH.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-NOTCH-SIGNALING")
dev.off()

pdf(paste0(name, "_escape_ridge_STAT3.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-IL6-JAK-STAT3-SIGNALING")
dev.off()


pdf(paste0(name, "_escape_ridge_MYOGENESIS.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-MYOGENESIS")
dev.off()


pdf(paste0(name, "_escape_ridge_APICAL-JUNCTION.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-APICAL-JUNCTION")
dev.off()

pdf(paste0(name, "_escape_ridge_APICAL-SURFACE.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-APICAL-SURFACE")
dev.off()


pdf(paste0(name, "_escape_ridge_TNFA.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-TNFA-SIGNALING-VIA-NFKB")
dev.off()

pdf(paste0(name, "_escape_ridge_E2F.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-E2F-TARGETS")
dev.off()


pdf(paste0(name, "_escape_ridge_ESTROGEN.pdf"))
ridgeEnrichment(scRep_example,
	        group.by="subclass",	
                assay = "escape.ssGSEA_normalized",
                gene.set = "HALLMARK-ESTROGEN-RESPONSE-LATE")
dev.off()



all.markers <- FindAllMarkers(scRep_example,
			      group.by="subclass",
                              assay = "escape.ssGSEA_normalized",
                              min.pct = 0,
                              logfc.threshold = 0)

df=data.frame(all.markers)
#saveRDS(scRep_example,file=paste0(name, "_escape.rds"))

write.table(df,file=paste0(name, "_escape_DEG.txt"),sep="\t",quote=F)

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
DefaultAssay(scRep_example)="escape.ssGSEA_normalized"

#pdf(paste0(name, "_escape_featureplot_angio.pdf"))
#FeaturePlot(scRep_example, "HALLMARK-ANGIOGENESIS") +
#  scale_color_gradientn(colors = colorblind_vector) +
#  theme(plot.title = element_blank())
#dev.off()

pdf(paste0(name, "_escape_featureplot_G2M.pdf"))
FeaturePlot(scRep_example, "HALLMARK-G2M-CHECKPOINT") +
  scale_color_gradient(colors = colorblind_vector) +
  theme(plot.title = element_blank())
dev.off()

pdf(paste0(name, "_escape_year.pdf"))
RidgePlot(
  scRep_example,
  "age_year",
  group.by = "subclass",
)
dev.off()
###, assay="escape.ssGSEA"
