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
@click.option('-l', '--label', type=click.STRING, default='majorclass', show_default=True, help='Label.')
@click.option('-W', '--width', type=click.FLOAT, default=16, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=5, show_default=True, help='Height.')
@click.option('-n', '--namefile', type=click.Path(exists=False, resolve_path=True), help='Feature name file.')
@click.option('-f', '--format', type=click.Choice(['pdf', 'svg']), is_flag=False, default='pdf', show_default=True, help='Image format.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, label, width, height, namefile, format, infile):
	"""
To plot gene-target enrichment heatmap.

INFILE: a .h5mu file.

\b
Example:
  f=scplusmdata.h5mu
  bname=ON_regulon
  outdir=$(mrrdir.sh)
  scenicplush5mu2eregulonheatmap -e scenicplus -d "$outdir" -b "$bname" -l majorclass -W 16 -H 5 -- "$f"

\b
Note:
  1. The variable "scRNA_counts:{label}" will be used for cell types.
  2. Feature name file ("-n|--namefile") is not headed with only one column for feature name subset.
  3. The resulting PDF can be editted by Acrobat but not Illustrator.

\b
See also:
  Upstream:
    ataccistopic2runscenicpluspipeline
  Related:
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
		f"label='{label}'",
		f"width={width}",
		f"height={height}",
		f"namefile='{namefile}'" if namefile else f"namefile=None",
		f"format='{format}'",
		f"infile='{infile}'",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
