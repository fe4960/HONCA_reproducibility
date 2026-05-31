#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exePython.exePython import exePython
from pathlib import Path
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-c', '--indexcol', type=click.INT, default=0, show_default=True, help='Index column.')
@click.option('-n', '--ntop', type=click.INT, default=10, show_default=True, help='Top regulon.')
@click.option('-m', '--numcolumn', type=click.INT, default=1, show_default=True, help='Number of columns in plot.')
@click.option('-W', '--width', type=click.FLOAT, default=5, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=5, show_default=True, help='Height.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, indexcol, ntop, numcolumn, width, height, infile):
	"""
To calculate regulon specificity scores using gene-based regulon enrichment scores.

INFILE: a .h5mu file.

\b
Example:
  f=$(mrrdir.sh ../scenicplush5mu2eregulonrss)/scplusmdata.txt.gz
  bname=ON_regulon
  outdir=$(mrrdir.sh)
  scenicplusrsstxt2plotrss -e scenicplus -d "$outdir" -b "$bname" -n 10 -- "$f"

\b
Note:
  1. Two modalities, ['direct_gene_based_AUC', 'extended_gene_based_AUC'], will be used calculation.

\b
See also:
  Upstream:
    ataccistopic2runscenicpluspipeline
    scenicplush5mu2eregulonrss
  Depends:
    Python/Scanpy

\b
Date: 2025/07/23
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f"{absdir}/python/{scriptname}.py"
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"indexcol={indexcol}",
		f"ntop={ntop}",
		f"numcolumn={numcolumn}",
		f"width={width}",
		f"height={height}",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
