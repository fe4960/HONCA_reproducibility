# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(jlutils))

pdf(NULL) # skip Rplots.pdf
addArchRThreads(threads=numthreads)

x=loadArchRProject(archrprjdir)
print('==> loadArchRProject()')
str(x)

markerfeature=getMarkerFeatures(
	ArchRProj=x,
	groupBy=group,
	useMatrix=matrixname,
	bias=c('TSSEnrichment', 'log10(nFrags)'),
	testMethod=testmethod,
	binarize=binarize,
	verbose=T
	)
print('==> getMarkerFeatures()')
str(markerfeature)
saveRDS(markerfeature, file=sprintf('%s/%s.rds', outdir, bname))
