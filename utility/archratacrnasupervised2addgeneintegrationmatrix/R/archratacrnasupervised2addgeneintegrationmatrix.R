# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(jlutils))

# RNA
if (endsWith(rnafile, '.rds')) {
	rna=readRDS(rnafile)
} else if (endsWith(rnafile, '.h5ad')) {
	rna=h5ad2seurat(rnafile, catobs=rnalabel)
} else {
	write(sprintf("Error: wrong format %s. Please input either .rds or .h5ad.\n", rnafile), stderr())
	q(status=1)
}
print('===> rna')
str(rna)

# ATAC
pdf(NULL) # skip Rplots.pdf
addArchRThreads(threads=numthreads)
atac=loadArchRProject(archrprjdir)
print('===> loadArchRProject()')
str(atac)

# Add integration between ATAC and RNA
atac=addGeneIntegrationMatrix(
	ArchRProj=atac,
	useMatrix=usematrix,
	matrixName=matrixname,
	reducedDims=dimreductname,
	seRNA=rna,
	groupATAC=ataclabel,
	groupRNA=rnalabel,
	groupList=NULL,
	sampleCellsATAC=10000,
	sampleCellsRNA=10000,
	nGenes=nhvg,
	reduction=anchorreduct,
	nameCell='predictedCell',
	nameGroup='predictedGroup',
	nameScore='predictedScore',
	addToArrow=T,
	force=T,
	verbose=T
)
print('===> addGeneIntegrationMatrix()')
str(atac)

saveArchRProject(
	ArchRProj=atac,
	outputDirectory=bname,
	load=F,
	overwrite=F
	)
write.txt(data.frame(cellNames=getCellNames(atac), atac@cellColData), file=sprintf('%s_metadata.txt.gz', bname))
