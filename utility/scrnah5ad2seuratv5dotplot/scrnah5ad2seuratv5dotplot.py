#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import click
import datetime
import shutil

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-W', '--width', type=click.FLOAT, default=8, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=8, show_default=True, help='Height.')
@click.option('-x', '--xlab', type=click.STRING, help='X-axis label.')
@click.option('-y', '--ylab', type=click.STRING, help='Y-axis label.')
@click.option('-g', '--groupby', type=click.STRING, default='leiden', show_default=True, help='Metadata to group cells.')
@click.option('-m', '--marker', type=click.Path(exists=True), required=True, help='A file of marker genes for cell types.')
@click.option('-c', '--ordergroup', type=click.Path(exists=True), help='The order of the group labels.')
@click.option('-n', '--norm', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Normalize raw counts.')
@click.option('-s', '--scale', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Scale and center the data.')
@click.option('-f', '--flip', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Flip axis.')
@click.option('--scalemin', type=click.FLOAT, help='Lower limit for scaling.')
@click.option('--scalemax', type=click.FLOAT, help='Upper limit for scaling.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-k', '--clean', is_flag=True, help='Clean temporary files.')
@click.option('-r', '--color', type=click.Choice(['red', 'blue']), is_flag=False, default='red', show_default=True, help='Color map.')
@click.argument('infile', type=click.Path(exists=True))
def main(outdir, bname, width, height, xlab, ylab, groupby, marker, ordergroup, norm, scale, flip, scalemin, scalemax, condaenv, clean, color, infile):
	"""
Generate a dotplot using SeuratV5::DotPlot() from a .h5ad file. This is to utilize BPCells and Seurat V5 to scale-upcalculation for an atlas-level .h5ad file. See scrnah5ad2seuratdotplot for conventional use case.

INFILE is a Scanpy object in a .h5ad format.

\b
Example:
  f=$(mrrdir.sh ../../../fixlabel/scrnah5ad2normscale/HRCA_AC.h5ad)
  labelmarker=../label_marker.txt
  labeltype=../label_type.txt
  group=AC_type
  bname=$(basename "$f" .h5ad)
  outdir=$(mrrdir.sh)
  scrnah5ad2seuratv5dotplot -e seurat -d "$outdir" -b "$bname" -m "$labelmarker" -g "$group" -c "$labeltype" -f -H 4 -W 24 -- "$f"

\b
Note:
  1. See "scrnah5ad2seuratdotplot -h" for usage.
  2. TODO: the timestamp of .h5ad file will change after BPCells::open_matrix_anndata_hdf5().
  3. TODO: BPCells functions are used to load sparse matrix files only, while metadata should be loaded separately.
  4. Timestamp bug is fixed in "jlutils::h5ad2seurat()" by adding flags="H5F_ACC_RDONLY" to "rhdf5::H5Fopen()".
  5. Seurat v5 operations are not broaded supported for on-disk data.
  https://github.com/satijalab/seurat/issues/8963#issuecomment-2142974481

\b
  ```
  Hi, Sorry, currently SCTransform() does not support on-disk data. For on-disk
  data (based on BPCells), we only support standard NormalizeData() function (as
  documented here:
  https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette).
  Improving SCTransform() for on-disk data would be our long-term goal.
  ```

\b
  6. Bug fix "-n|--norm".

\b
  ```Error message
  > Matrix(LayerData(x, assay='RNA', layer='data'), sparse=T)
  Error in as.vector(data) :
    no method for coercing this S4 class to a vector
  ```

\b
  ```Bug fix
  - counts=Matrix(LayerData(x, assay='RNA', layer='data'), sparse=T)
  + counts=LayerData(x, assay='RNA', layer='data')
    ```

\b
See also:
  Related:
    scrnah5ad2seuratdotplot
    scrnah5adnsforestmarkertxt2seuratv5dotplot
  Depends:
    R/SeuratV5
    R/BPCells

\b
Date: 2024/06/04
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)

	# Copy infile .h5ad because timestamp will be changed after BPCells::open_matrix_anndata_hdf5(). TODO
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')
	tmpdir=f"{outdir}/tmpdir_{bname}_{nowtimestr}"
	tmpfile=f"{tmpdir}/{bname}.h5ad"
	bpcdir=f"{tmpdir}/BPCells_{bname}"
	Path(tmpdir).mkdir(parents=True, exist_ok=True)
	shutil.copy2(infile, tmpfile)

	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/R/{scriptname}.R"
	exprs=[
		f"infile='{tmpfile}'",
		f"bpcdir='{bpcdir}'",
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"width={width}",
		f"height={height}",
		f"xlab='{xlab}'" if xlab else f"xlab=NULL",
		f"ylab='{ylab}'" if ylab else f"ylab=NULL",
		f"groupby='{groupby}'",
		f"markerfile='{marker}'",
		f"ordergroup='{ordergroup if ordergroup is not None else ''}'",
		f"norm={norm}",
		f"scale={scale}",
		f"flip={flip}",
		f"scalemin={scalemin}" if scalemin else f"scalemin=NA",
		f"scalemax={scalemax}" if scalemax else f"scalemax=NA",
		f"color='{color}'",
		]
	exitstatus=exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

	if clean:
		# clear tmpdir
		shutil.rmtree(tmpdir)

	return exitstatus

if __name__ == "__main__":
	sys.exit(main())
