# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(RColorBrewer))

pdf(NULL) # skip Rplots.pdf
addArchRThreads(threads=numthreads)
atac=loadArchRProject(archrprjdir)
print('===> loadArchRProject()')
str(atac)

if (!is.null(key) && length(value)>0) {
	if (invert) {
		cellsKeep=atac$cellNames[!(atac@cellColData[[key]] %in% value)]
	} else {
		cellsKeep=atac$cellNames[atac@cellColData[[key]] %in% value]
	}

	atac=subsetArchRProject(
		ArchRProj=atac,
		cells=cellsKeep,
		outputDirectory=paste0(bname, "_subset"),
		dropCells=T,
		)

	print('===> subset()')
	str(atac)
}

# set color
if (setcolor) {
	metagroup=as.data.frame(atac@cellColData)[[group]]

	if (is.factor(metagroup) || is.character(metagroup)) {
		grouplabel=if (is.factor(metagroup)) { levels(metagroup) } else { unique(metagroup) }
		colornames=c(
			brewer.pal(9, "Set1"),
			brewer.pal(8, "Set2"),
			brewer.pal(12, "Set3")
			)
		color_map=setNames(colornames[seq_along(grouplabel)], sort(grouplabel))
	} else {
		color_map=NULL
	}
	str(color_map)
	print(color_map)
}

p=plotPeak2GeneHeatmap(
	ArchRProj=atac,
	groupBy=group,
	corCutOff=corcutoff,
	FDRCutOff=fdrcutoff,
	varCutOffATAC=0.25,
	varCutOffRNA=0.25,
	k=knn,
	nPlot=25000,
	palGroup=if (setcolor) color_map else NULL,
	verbose=T
	)

if (ggplot) {
	pdf(sprintf('%s.pdf', bname), width=width, height=height)
	ComplexHeatmap::draw(p)
	dev.off()
} else {
	plotPDF(
		p,
		name=sprintf('%s.pdf', bname),
		width=width,
		height=height,
		addDOC=F
		)
}
