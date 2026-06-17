library(ArchR)

dir=paste0("/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/ON_ATAC/")
cell="major_full"
proj2=loadArchRProject(path = dir, force = FALSE, showLogo = TRUE)

peakset=getPeakSet(proj2)
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/"

saveRDS(peakset,file=paste0(dir1,"/peakset_clusters_",cell,".rds"))
