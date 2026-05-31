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
@click.option('-u', '--usematrix', type=click.Choice(['GeneScoreMatrix', 'TileMatrix']), is_flag=False, default='GeneScoreMatrix', show_default=True, help='Matrix for integration.')
@click.option('-m', '--matrixname', type=click.STRING, default='GeneIntegrationMatrix', show_default=True, help='Matrix name.')
@click.option('-c', '--dimreductname', type=click.Choice(['IterativeLSI', 'Harmony']), is_flag=False, default='Harmony', show_default=True, help='Reduced dimension.')
@click.option('-r', '--rnafile', type=click.Path(exists=True, resolve_path=True), help='RNA file.')
@click.option('-l', '--rnalabel', type=click.STRING, default='majorclass', show_default=True, help='RNA label.')
@click.option('-L', '--ataclabel', type=click.STRING, default='majorclass', show_default=True, help='ATAC label.')
@click.option('-n', '--nhvg', type=click.INT, default=10000, show_default=True, help='Number of highly variable genes.')
@click.option('-a', '--anchorreduct', type=click.Choice(['cca', 'rpca', 'lsiproject', 'pcaproject']), is_flag=False, default='cca', show_default=True, help='Dimension reduction to find anchors.')
@click.argument('archrprjdir', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, numthreads, usematrix, matrixname, dimreductname, rnafile, rnalabel, ataclabel, nhvg, anchorreduct, archrprjdir):
	"""
To integrate an ArchR project with a Seurat object using constrained/superviserd integration.

ARCHRPRJDIR: an ArchR project directory.

\b
Example:
  atac=$(mrrdir.sh ../archrproject2addpeakset)/ON_ATAC
  rna=/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/ONONH/coembed/dotplot/woRPE/scrnah5adaddmetadatamerge/ON_RNA.h5ad
  outdir=$(mrrdir.sh)
  function cmd {
  local nhvg=$1
  local bname=ON_ATAC_RNA_${nhvg}
  slurmucihpc.sh -b -t 2 -m 400G -w 3-00:00:00 -- archratacrnasupervised2addgeneintegrationmatrix -e r-archr -d "$outdir" -b "$bname" -r "$rna" -l majorclass -L majorclass -c Harmony -n "$nhvg" -a cca -- "$atac"
  }
  source env_parallel.bash
  env_parallel -j 1 cmd ::: 2000 10000

\b
Note:
  1. To perform ArchR analysis on scATAC-Seq data.
  https://www.archrproject.com/articles/Articles/tutorial.html
  https://www.archrproject.com/bookdown/cross-platform-linkage-of-scatac-seq-cells-with-scrna-seq-cells.html
  https://www.archrproject.com//bookdown/peak2genelinkage-with-archr.html
  2. The "-a|--anchorreduct" option can only use "cca". Other dimension reductions have not yet been implemented in ArchR.
  3. archratacrnaintegrate2peak2genelinks is split into two wrappers, archratacrnasupervised2addgeneintegrationmatrix and archrproject2addpeak2genelinks.
  3.1 The integration step is too slow in archratacrnaintegrate2peak2genelinks.
  4. Intermediate RNAintegration/ is stored in the input ArchR object folder. Therefore, do not run multiple instances of this wrapper in parallel.

\b
See also:
  Related:
    archratacrnaintegrate2peak2genelinks
    archrseuratintegrateunconstrain2addgeneintegrationmatrix
  Downstream:
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
		f"usematrix='{usematrix}'",
		f"matrixname='{matrixname}'",
		f"dimreductname='{dimreductname}'",
		f"rnafile='{rnafile}'",
		f"rnalabel='{rnalabel}'",
		f"ataclabel='{ataclabel}'",
		f"nhvg={nhvg}",
		f"anchorreduct='{anchorreduct}'",
		f"archrprjdir='{archrprjdir}'",
		]
	os.chdir(outdir)
	return exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
