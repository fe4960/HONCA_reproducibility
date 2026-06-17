library(escape)
library(Seurat)
library(BiocParallel)
library(GSVA)
library(msigdbr)
library(dplyr)
#name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat"
args <- commandArgs(trailingOnly = TRUE)


#Sys.setenv(PATH = paste(
#  "/pub/junw42/local/miniconda/envs/escape/bin",
#  Sys.getenv("PATH"),
#  sep = ":"
#))

Sys.setenv(R_SCRIPT = "/pub/junw42/local/miniconda/envs/escape/bin/Rscript")

name=args[1]

scRep_example=readRDS(paste0(name, ".rds"))
path=paste0(name, ".obs.gz")
if(file.exists(path)){
meta=read.table(path, header=T, sep=",", row.names=1)
scRep_example <- AddMetaData(scRep_example, metadata = meta)
}

#scRep_example=readRDS(paste0(name, "_simple.rds"))
#GS.hallmark <- getGeneSets(library = "H", subcategory="GO:BP")
######GS.hallmark <- getGeneSets(library = "H",species="Homo sapiens") #, subcategory="GO:BP") #%>% split(x = .$gene_symbol, f = .$gs_name)
NSC_quiescent <- c(
  "GFAP", "SOX2", "NES", "VIM", "FABP7",
  "HES1", "HES5", "ID1", "ID2", "ID3",
  "PROM1", "AQP4", "SLC1A3", "ALDH1L1", "GLI3"
)

NSC_activated <- c(
  "SOX2", "NES", "VIM", "EGFR", "ASCL1",
  "MKI67", "TOP2A", "PCNA", "MCM2", "MCM6",
  "HMGB2", "TUBB", "HES1", "DLL1", "CCND1"
)

NPC_early <- c(
  "ASCL1", "DLL1", "EGFR", "SOX4", "SOX11",
  "DCX", "ELAVL4", "STMN2", "TUBB3", "NEUROD1",
  "NEUROD2", "HES6", "BTG2", "MKI67", "TOP2A"
)

Astrocyte_homeostatic <- c(
  "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1",
  "FGFR3", "ATP1A2", "GLUL", "APOE", "SOX9"
)

GS_NSC <- list(
  NSC_quiescent = NSC_quiescent,
  NSC_activated = NSC_activated,
  NPC_early = NPC_early,
  Astrocyte_homeostatic = Astrocyte_homeostatic
)



param <- BiocParallel::SnowParam(workers = 1)

#Sys.setenv(PATH = paste("/pub/junw42/local/miniconda/envs/escape/bin", Sys.getenv("PATH"), sep=":"))
#stopifnot(nzchar(Sys.which("Rscript")))
#system("Rscript --version")
#Sys.which("R")
#Sys.which("Rscript")

#cat("R      :", Sys.which("R"), "\n")
#cat("Rscript:", Sys.which("Rscript"), "\n")
#cat("R_HOME :", Sys.getenv("R_HOME"), "\n")
#cat("R.home :", R.home(), "\n")
#cat("TMPDIR :", Sys.getenv("TMPDIR"), "\n")
#cat("tempdir:", tempdir(), "\n")

system2(Sys.which("Rscript"), "--version")

scRep_example <- runEscape(input.data=scRep_example, 
                           gene.sets = GS_NSC, #GS.hallmark,
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
                                      gene.sets = GS_NSC ) #GS.hallmark)
#				      msigdbGeneSets = c("H"),
#				      )

scRep_example <- performNormalization(scRep_example,
                                      assay = "escape.ssGSEA",
		#		      method= "ssGSEA",
		#		      min.size=5,
                                     gene.sets = GS_NSC, #GS.hallmark,
#				      msigdbGeneSets = c("H"),
                                      scale.factor = scRep_example$nFeature_RNA)

#scRep_example$nFeature_RNA

saveRDS(scRep_example,file=paste0(name, "_escape_NSC.rds"))

library(ggplot2)
library(RColorBrewer)
pdf(paste0(name, "_escape_heatmap_NSC.pdf"), width=12, height = 8)
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = "all", #rownames(scRep_example@assays$escape.ssGSEA@layers$data)[1:15],
		  #                  gene.set.use = rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.ssGSEA",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))+scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 
dev.off()


pdf(paste0(name, "_escape_heatmap_norm_NSC.pdf"), width=12, height = 8)
heatmapEnrichment(scRep_example,
                  group.by = "subclass",
                  gene.set.use = "all", #rownames(scRep_example@assays$escape.ssGSEA_normalized@layers$data)[1:15],
		  #                  gene.set.use = rownames(scRep_example@assays$escape.ssGSEA@data)[1:15],
                  assay = "escape.ssGSEA_normalized",
		  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))+scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) 
dev.off()





all.markers <- FindAllMarkers(scRep_example,
			      group.by="subclass",
                              assay = "escape.ssGSEA_normalized",
                              min.pct = 0,
                              logfc.threshold = 0)

df=data.frame(all.markers)

write.table(df,file=paste0(name, "_escape_DEG_NSC.txt"),sep="\t",quote=F)

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
DefaultAssay(scRep_example)="escape.ssGSEA_normalized"


#pdf(paste0(name, "_escape_featureplot_G2M.pdf"))
#FeaturePlot(scRep_example, "HALLMARK-G2M-CHECKPOINT") +
#  scale_color_gradientn(colors = colorblind_vector) +
#  theme(plot.title = element_blank())
#dev.off()

pdf(paste0(name, "_escape_year_NSC.pdf"))
RidgePlot(
  scRep_example,
  "age_year",
  group.by = "subclass",
)
dev.off()
