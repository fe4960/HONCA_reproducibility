#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import random
import socket
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-0', '--empty', type=click.Choice(['F', 'T']), is_flag=False, flag_value='F', default='T', show_default=True, help='Use empty count matrix.')
@click.option('-g', '--group', type=click.STRING, multiple=True, required=True, help='Group to plot the UMAP.')
@click.option('-n', '--numericgroup', type=click.STRING, multiple=True, help='Numeric group to plot the UMAP.')
@click.option('-s', '--split', type=click.STRING, default='', show_default=True, help='Split name.')
@click.option('-H', '--height', type=click.FLOAT, default=6, show_default=True, help='Height.')
@click.option('-W', '--width', type=click.FLOAT, default=6, show_default=True, help='Width.')
@click.option('-L', '--legendwidth', type=click.FLOAT, default=1, show_default=True, help='Legend width.')
@click.option('-t', '--showtitle', type=click.Choice(['F', 'T']), is_flag=False, flag_value='F', default='T', show_default=True, help='Keep plot title.')
@click.option('-a', '--noaxis', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='No axis.')
@click.option('-f', '--format', type=click.Choice(['pdf', 'png']), is_flag=False, flag_value='pdf', default='png', show_default=True, help='Image format.')
@click.option('-c', '--setcolor', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Set colors.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, empty, group, numericgroup, split, height, width, legendwidth, showtitle, noaxis, format, setcolor, infile):
	"""
Plot UMAP by Seurat::DimPlot().

INFILE is a .h5ad file with UMAP.

\b
Example:
  indir=$(mrrdir.sh ../glueh5adcoembed2labeltransfer)
  outdir=$(mrrdir.sh)
  function cmd {
  local f=$1
  local bname=$(basename "$f" .h5ad)
  if fileexists.sh "$f"
  then
  	slurmucihpc.sh -p free-gpu -g 0 -t 2 -m 40G -w 2:00:00 -- scrnah5ad2seuratumapby -e seuratv4 -d "$outdir" -b "$bname" -s scglue_domain__ -g scglue_celltype__ -g scglue_prediction -- "$f"
  fi
  }
  source env_parallel.bash
  env_parallel cmd ::: "$indir"/*.h5ad

\b
Note:
  1. For "-0|--empty", only .obs and .obsm are loaded from the .h5ad file. This is useful for .h5ad files without .X or for efficient loading. This feature has been tested with Seurat v4. TODO: Test with Seurat v5.

\b
See also:
  Related:
    saturnh5adumap2seuratdimplotsplitby
    scrnaseuratumapby.sh
    scrnah5adumapby.sh
  Depends:
    R/jlutils
    R/Seurat

\b
Date: 2025/01/10
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"empty={empty}",
		f"group=c({str(list(group))[1:-1] if group is not None else ''})",
		f"numericgroup=c({str(list(numericgroup))[1:-1] if numericgroup is not None else ''})",
		f"split='{split}'",
		f"height={height}",
		f"width={width}",
		f"legendwidth={legendwidth}",
		f"showtitle={showtitle}",
		f"noaxis={noaxis}",
		f"format='{format}'",
		f"setcolor={setcolor}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
