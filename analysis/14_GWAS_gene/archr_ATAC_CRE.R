library(ArchR)
set.seed(1)
addArchRThreads(threads = 10)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
addArchRGenome("hg38")
#cell=args[1]
#cell="final_major_full_max"
#dir=paste0("/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_",cell)

#dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/ON_ATAC/")
dir=paste0("/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/ON_ATAC/")
proj5=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)


proj5 <- addPeak2GeneLinks(
    ArchRProj = proj5,
    reducedDims = "IterativeLSI",
    maxDist = 250000
)

