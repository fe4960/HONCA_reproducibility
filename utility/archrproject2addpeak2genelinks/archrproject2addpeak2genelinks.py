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
@click.option('-m', '--matrixname', type=click.STRING, default='GeneIntegrationMatrix', show_default=True, help='Matrix name.')
@click.option('-c', '--dimreductname', type=click.Choice(['IterativeLSI', 'Harmony']), is_flag=False, default='Harmony', show_default=True, help='Reduced dimension.')
@click.argument('archrprjdir', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, numthreads, matrixname, dimreductname, archrprjdir):
	"""
To add peak2gene links from ATAC-RNA integration matrix from archratacrnasupervised2addgeneintegrationmatrix.

ARCHRPRJDIR: an ArchR project directory.

\b
Example:
  indir=$(mrrdir.sh ..)
  outdir=$(mrrdir.sh)
  function cmd {
  local bname=$1
  local archrdir=$indir/$bname
  slurmucihpc.sh -t 2 -m 100G -w 3-00:00:00 -- archrproject2addpeak2genelinks -e r-archr -d "$outdir" -b "$bname" -m GeneIntegrationMatrix -c Harmony -- "$archrdir"
  }
  source env_parallel.bash
  env_parallel cmd ::: ON_ATAC_RNA_2000 ON_ATAC_RNA_10000

\b
Note:
  1. To perform ArchR analysis on scATAC-Seq data.
  https://www.archrproject.com/articles/Articles/tutorial.html
  https://www.archrproject.com/bookdown/cross-platform-linkage-of-scatac-seq-cells-with-scrna-seq-cells.html
  https://www.archrproject.com//bookdown/peak2genelinkage-with-archr.html

\b
See also:
  Related:
    archratacrnaintegrate2peak2genelinks
  Upstream:
    archratacrnasupervised2addgeneintegrationmatrix
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
		f"matrixname='{matrixname}'",
		f"dimreductname='{dimreductname}'",
		f"archrprjdir='{archrprjdir}'",
		]
	os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
