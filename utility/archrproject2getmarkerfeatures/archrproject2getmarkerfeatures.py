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
@click.option('-g', '--group', type=click.STRING, default='ClustersHarmony', show_default=True, help='Group.')
@click.option('-n', '--matrixname', type=click.Choice(['GeneScoreMatrix', 'TileMatrix', 'PeakMatrix']), is_flag=False, default='GeneScoreMatrix', show_default=True, help='Matrix name.')
@click.option('-m', '--testmethod', type=click.Choice(['wilcoxon', 'ttest', 'binomial']), is_flag=False, default='wilcoxon', show_default=True, help='Test name.')
@click.option('-z', '--binarize', type=click.Choice(['F', 'T']), is_flag=False, flag_value='T', default='F', show_default=True, help='Binarize the matrix.')
@click.argument('archrprjdir', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, numthreads, group, matrixname, testmethod, binarize, archrprjdir):
	"""
To identify marker features from an ArchR project.

ARCHRPRJDIR: an ArchR project directory.

\b
Example:
  archrprjdir=$(mrrdir.sh ../archrcellcoldata2metadatabybarcode)/archr
  bname=archr
  outdir=$(mrrdir.sh)
  slurmucihpc.sh -t 4 -m 80G -w 20:00:00 -p free-gpu -g 0 -- archrproject2getmarkerfeatures -e r-archr -d "$outdir" -b "$bname" -g ClustersHarmony -n GeneScoreMatrix -m wilcoxon -- "$archrprjdir"
  # @Deprecated. Get markers. See archrmarkerfeaturerds2marker
  archrproject2getmarkerfeatures -e r-archr -d "$outdir" -b "$bname" -g ClustersHarmony -n GeneScoreMatrix -m wilcoxon -N 20 -- "$archrprjdir"
  archrproject2getmarkerfeatures -e r-archr -d "$outdir" -b "$bname" -g ClustersHarmony -n GeneScoreMatrix -m wilcoxon -c "FDR <= 0.01 & Log2FC >= 0.5" -- "$archrprjdir"

\b
Note:
  1. To perform ArchR analysis on scATAC-Seq data.
  https://www.archrproject.com/articles/Articles/tutorial.html
  2. This wrapper can be used to call DARs using the PeakMatrix.
  3. Deprecated get markers. See archrmarkerfeaturerds2marker.

\b
  ```
  -@click.option('-c', '--cutoff', type=click.STRING, help='Cut off.')
  -@click.option('-N', '--ntop', type=click.INT, show_default=True, help='Top ranked genes.')
  ```

\b
See also:
  Downstream:
    archrmarkerfeaturerds2marker
    archrmarkerfeature2heatmap
  Depends:
    R/ArchR

\b
Date: 2025/05/24
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
		f"matrixname='{matrixname}'",
		f"testmethod='{testmethod}'",
		f"binarize={binarize}",
		f"archrprjdir='{archrprjdir}'",
		]
	os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
