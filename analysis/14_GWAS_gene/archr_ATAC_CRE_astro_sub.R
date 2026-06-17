library(ArchR)
#library(Seurat)
library(GenomicRanges)
set.seed(1)
addArchRThreads(threads = 20)
library(tidyr)
library(parallel)


#dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/7_atac/project_rnav2_downsample10k_latent30_dim30/")
#dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/7_atac/project_rnav2_downsample10k_latent30_dim30_majortype_transfer/")

#dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/ON_ATAC/")
#dir=paste0("/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/archr_object/ON_ATAC_RNA_10000/")
dir=paste0("/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/astrocyte_archr/type4/ON_astrocyte/")

#peak_r=GRanges(seqnames="chr3",ranges=IRanges(56337404, 57337404))

#peak_r=GRanges(seqnames="chr3",ranges=IRanges(56537404, 57137404))

#peak_r=GRanges(seqnames="chr9",ranges=IRanges(21749657, 22349657))
peak_r=GRanges(seqnames="chr9",ranges=IRanges(21849657, 22249657))
#    name = "Plot-Tracks-DAR_g2p_CDKN2B_majorclass.pdf",


proj5=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

p <- plotBrowserTrack(
    ArchRProj = proj5,
    region=peak_r,
    #useGroups=cells,
    groupBy = "celltype",
 loops = getPeak2GeneLinks(proj5,resolution=1,corCutOff = 0.6,FDRCutOff = 0.01)
#  loops = getCoAccessibility(proj5,corCutOff = 0.5,resolution = 1,returnLoops = TRUE)
######## loops = getPeak2GeneLinks(proj5,resolution=1,corCutOff = 0.6,FDRCutOff = 0.01)
)


plotPDF( p, ################plotList = p,
#   name = "Plot-Tracks-DAR_g2p_ARHGEF3_astro_sub.pdf",
#    name= "Plot-Tracks-DAR_g2p_CDKN2B_astro_sub.pdf",
    name= "Plot-Tracks-DAR_g2p_CDKN2B_astro_sub_200k.pdf",

    ArchRProj = proj5,
    addDOC = FALSE, width = 8, height = 7)

