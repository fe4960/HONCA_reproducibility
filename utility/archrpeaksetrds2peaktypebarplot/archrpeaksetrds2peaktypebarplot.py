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
@click.option('-W', '--width', type=click.FLOAT, default=3.5, show_default=True, help='Width.')
@click.option('-H', '--height', type=click.FLOAT, default=3.5, show_default=True, help='Height.')
@click.option('-p', '--position', type=click.Choice(['stack', 'dodge']), is_flag=False, flag_value='dodge', default='stack', show_default=True, help='Position.')
@click.option('-l', '--legendposition', type=click.Choice(['left', 'top', 'right', 'bottom']), is_flag=False, flag_value='bottom', default='right', show_default=True, help='Legend position.')
@click.option('-x', '--xlab', type=click.STRING, help='Group.')
@click.option('-y', '--ylab', type=click.STRING, help='Percentage.')
@click.option('-m', '--main', type=click.STRING, help='Main.')
@click.option('-s', '--size', type=click.INT, default=5, show_default=True, help='Axis text size.')
@click.option('-M', '--mode', type=click.Choice(['proportion', 'number']), is_flag=False, flag_value='proportion', default='number', show_default=True, help='Number or proportion.')
@click.option('-g', '--group', type=click.STRING, multiple=True, help='Include only these Group labels (repeatable). By default all groups are included.')
@click.option('-n', '--invert', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Invert group selection.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, width, height, position, legendposition, xlab, ylab, main, size, mode, group, invert, infile):
	"""
To plot the peak summary from the calculated peak set using archrproject2addpeakset.

INFILE: the peakset .rds file.

\b
Example:
  infile=/tmp/archr/ON_ATAC_peakset.rds
  bname=$(basename "$infile" _peakset.rds)
  outdir=$(mktemp -d -u)
  archrpeaksetrds2peaktypebarplot -d "$outdir" -b "$bname" -y 'Number of OCRs' -- "$infile"
  archrpeaksetrds2peaktypebarplot -d "$outdir" -b "$bname" -y 'Proportion' -M -- "$infile"
  archrpeaksetrds2peaktypebarplot -d "$outdir" -b "$bname" -g 'GroupA' -g 'GroupB' -- "$infile"
  archrpeaksetrds2peaktypebarplot -d "$outdir" -b "$bname" -g 'MG(n = 17139)' -n T -- "$infile"

\b
Note:
  1. To perform ArchR analysis on scATAC-Seq data.
  https://www.archrproject.com//bookdown/calling-peaks-w-macs2.html

\b
See also:
  Upstream:
    archrproject2addpeakset

\b
Date: 2025/05/22
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"width={width}",
		f"height={height}",
		f"position='{position}'",
		f"legendposition='{legendposition}'",
		f"xlab='{xlab if xlab is not None else ''}'",
		f"ylab='{ylab if ylab is not None else ''}'",
		f"main='{main if main is not None else ''}'",
		f"size={size}",
		f"mode='{mode}'",
		f"group=c({','.join(repr(g) for g in group) if group else ''})",
		f"invert={invert}",
		f"infile='{infile}'",
		]
	os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
