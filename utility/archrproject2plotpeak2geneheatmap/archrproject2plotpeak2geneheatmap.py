#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-t', '--numthreads', type=click.INT, default=12, show_default=True, help='Number of threads.')
@click.option('-g', '--group', type=click.STRING, default='majorclass', show_default=True, help='Group.')
@click.option('-c', '--corcutoff', type=click.FLOAT, default=0.5, show_default=True, help='Min peak-to-gene correlation.')
@click.option('-f', '--fdrcutoff', type=click.FLOAT, default=0.01, show_default=True, help='Max peak-to-gene FDR.')
@click.option('-k', '--knn', type=click.INT, default=30, show_default=True, help='KNN to group peak-to-gene links.')
@click.option('-W', '--width', type=click.FLOAT, default=10, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=10, show_default=True, help='Height.')
@click.option('-G', '--ggplot', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Save by ggplot2::ggsave() or ArchR::plotPDF().')
@click.option('-s', '--setcolor', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Set colors.')
@click.option('-K', '--key', type=click.STRING, help='Key to subset by (not working) (TODO).')
@click.option('-V', '--value', type=click.STRING, multiple=True, help='Values for the subset key (repeatable).')
@click.option('-n', '--invert', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Invert group selection.')
@click.argument('archrprjdir', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, numthreads, group, corcutoff, fdrcutoff, knn, width, height, ggplot, setcolor, key, value, invert, archrprjdir):
	"""
To plot side by side heatmaps of linked ATAC and Gene regions from addPeak2GeneLinks().

ARCHRPRJDIR: ArchR object folder.

\b
Example:
  indir=$(mrrdir.sh ../archrproject2addpeak2genelinks)
  outdir=$(mrrdir.sh)
  function cmd {
  local archrname=$1
  local f=$indir/$archrname
  local ggsave=$2
  local knn=$3
  local bname=${archrname}_${knn}_${ggsave}
  slurmucihpc.sh -t 2 -m 80G -w 2:00:00 -- archrproject2plotpeak2geneheatmap -e r-archr -d "$outdir" -b "$bname" -g majorclass -c 0.5 -f 0.01 -k "$knn" -W 10 -H 10 -G "$ggsave" -- "$f"
  }
  source env_parallel.bash
  env_parallel cmd ::: ON_ATAC_RNA_2000 ON_ATAC_RNA_10000 ::: F T ::: 25 30 40 50 80 100

\b
Note:
  1. To perform ArchR analysis on scATAC-Seq data.
  https://www.archrproject.com/articles/Articles/tutorial.html

\b
See also:
  Upstream:
    archratacrnaintegrate2addpeak2genelinks
    archratacrnasupervised2addgeneintegrationmatrix
    archrproject2addpeak2genelinks
  Depends:
    R/ArchR

\b
Date: 2025/05/25
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"numthreads={numthreads}",
		f"group='{group}'",
		f"corcutoff={corcutoff}",
		f"fdrcutoff={fdrcutoff}",
		f"knn={knn}",
		f"width={width}",
		f"height={height}",
		f"ggplot={ggplot}",
		f"setcolor={setcolor}",
		f"key='{key}'" if key is not None else f"key=NULL",
		f"value=c({','.join(repr(v) for v in value) if value else ''})",
		f"invert={invert}",
		f"archrprjdir='{archrprjdir}'",
		]
	os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
