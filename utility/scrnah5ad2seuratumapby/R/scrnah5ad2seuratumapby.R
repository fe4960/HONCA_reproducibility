# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(Seurat))

if (empty) {
	indata=h5ad2seurat0(infile, catobs=c(group, split), numobs=numericgroup, obsm='X_umap')
} else {
	indata=h5ad2seurat(infile, catobs=c(group, split), numobs=numericgroup, obsm='X_umap')
}
print('==> indata')
str(indata)

# bug fix points for PDF
raster=NULL
if (format=='pdf') {
	raster=F # to save all vectorized points
}

if (setcolor) {
	colornames=c(
		brewer.pal(9, "Set1")
		, brewer.pal(8, "Set2")
		, brewer.pal(12, "Set3")
		)
} else {
	colornames=NULL
}

mclapply(
	c(group, if(nzchar(split)) split else NULL),
	function(g) {
		metagroup=indata@meta.data[[g]]
		color_map=if (setcolor && is.factor(metagroup) || is.character(metagroup)) {
			grouplabel=levels(metagroup) %||% unique(metagroup)
			setNames(colornames[seq_along(grouplabel)], sort(grouplabel))
		} else NULL
		str(color_map)
		print(color_map)

		seuratdimplotgroupby(indata, reduct='umap', group=g, label=F, outfile=sprintf('%s/%s_umap_%s_wolabel.%s', outdir, bname, g, format), width=width+legendwidth, height=height, nolegend=F, noaxis=noaxis, showtitle=showtitle, raster=raster, colors=color_map)
		seuratdimplotgroupby(indata, reduct='umap', group=g, label=T, outfile=sprintf('%s/%s_umap_%s_ondata.%s', outdir, bname, g, format), width=width, height=height, nolegend=T, noaxis=noaxis, showtitle=showtitle, raster=raster, colors=color_map)
	})

if (nzchar(split)) {
	nsplit=nlevels(indata@meta.data[[split]])
	str(nsplit)
	mclapply(
		group,
		function(g) {
			metagroup=indata@meta.data[[g]]
			color_map=if (setcolor && is.factor(metagroup) || is.character(metagroup)) {
				grouplabel=levels(metagroup) %||% unique(metagroup)
				setNames(colornames[seq_along(grouplabel)], sort(grouplabel))
			} else NULL
			str(color_map)
			print(color_map)

			seuratdimplotsplitby(indata, reduct='umap', group=g, split=split, label=F, outfile=sprintf('%s/%s_%s_%s_wolabel.%s', outdir, bname, split, g, format), width=width*nsplit+legendwidth, height=height, nolegend=F, noaxis=noaxis, showtitle=showtitle, raster=raster, colors=color_map)
			seuratdimplotsplitby(indata, reduct='umap', group=g, split=split, label=T, outfile=sprintf('%s/%s_%s_%s_ondata.%s', outdir, bname, split, g, format), width=width*nsplit, height=height, nolegend=T, noaxis=noaxis, showtitle=showtitle, raster=raster, colors=color_map)
		})
}
