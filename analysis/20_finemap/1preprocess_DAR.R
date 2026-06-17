library(ArchR)
#/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/DAR/ON_ATAC_binary_marker.rds
dir="/dfs3b/ruic20_lab/jinl14/share/HECA/ONONH/atac/DAR/"
cell="ON_ATAC"
markersPeaks=readRDS(paste0(dir,"/",cell,"_binary_marker.rds"))
#markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
#/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/DAR/
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/DAR/"
for(c in 1:length(markersPeaks)){
file=paste0(dir1,"/",names(markersPeaks[c]),"_DAR_peak_FDR005")
peak=data.frame(markersPeaks[[c]])
peak=peak[(peak$FDR <=0.05)&(peak$Log2FC >= 0.5),]
write.table(peak,file=file,quote=F,sep="\t")
}
#write.table(markerList[[c]],file=file,quote=F,sep="\t")
#}
