library(ArchR)
#library(Seurat)
library(GenomicRanges)
set.seed(1)
addArchRThreads(threads = 20)
library(tidyr)
library(parallel)


#dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/7_atac/project_rnav2_downsample10k_latent30_dim30/")
#dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/7_atac/project_rnav2_downsample10k_latent30_dim30_majortype_transfer/")

dir=paste0("/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/ON_ATAC/")

peak_r=GRanges(seqnames="chr9",ranges=IRanges(21749657, 22349657))


#peak_r=GRanges(seqnames="chr3",ranges=IRanges(56337404, 57337404))

proj5=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)
p <- plotBrowserTrack(
    ArchRProj = proj5,
    region=peak_r,
    #useGroups=cells,
    groupBy = "majorclass",
 loops = getPeak2GeneLinks(proj5,resolution=1,corCutOff = 0.5,FDRCutOff = 0.01)
#  loops = getCoAccessibility(proj5,corCutOff = 0.5,resolution = 1,returnLoops = TRUE)
######## loops = getPeak2GeneLinks(proj5,resolution=1,corCutOff = 0.6,FDRCutOff = 0.01)
)

plotPDF( p, ################plotList = p,
   name = "Plot-Tracks-DAR_g2p_CDKN2B_majorclass.pdf",

	#   name = "Plot-Tracks-DAR_g2p_ARHGEF3_majorclass.pdf",
    ArchRProj = proj5,
    addDOC = FALSE, width = 8, height = 7)

