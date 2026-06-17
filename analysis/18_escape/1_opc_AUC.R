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


pdf(paste0(name, "_escape_heatmap.pdf"))
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.GSVA",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE,
		  palette = "Spectral")
dev.off()


all.markers <- FindAllMarkers(scRep_example,
                              assay = "escape.ssGSEA_normalized",
                              min.pct = 0,
                              logfc.threshold = 0)

df=data.frame(all.markers)
saveRDS(scRep_example,file=paste0(name, "_escape.rds"))

write.table(df,file=paste0(name, "_escape_DEG.txt"),sep="\t",quote=F)
