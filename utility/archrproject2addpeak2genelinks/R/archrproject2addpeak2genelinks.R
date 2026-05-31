# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(jlutils))

# ATAC
pdf(NULL) # skip Rplots.pdf
addArchRThreads(threads=numthreads)
atac=loadArchRProject(archrprjdir)
print('===> loadArchRProject()')
str(atac)

# Add peak to gene links
atac=addPeak2GeneLinks(
	ArchRProj=atac,
	reducedDims=dimreductname,
	useMatrix=matrixname,
	maxDist=250000,
	verbose=T
)
print('===> addPeak2GeneLinks()')
str(atac)

saveArchRProject(
	ArchRProj=atac,
	outputDirectory=bname,
	load=F,
	overwrite=F
	)
write.txt(data.frame(cellNames=getCellNames(atac), atac@cellColData), file=sprintf('%s_metadata.txt.gz', bname))
