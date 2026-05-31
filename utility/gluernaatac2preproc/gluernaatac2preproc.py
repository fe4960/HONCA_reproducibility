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
@click.option('-r', '--rna', type=click.Path(exists=True, resolve_path=True), required=True, help='RNA.')
@click.option('-a', '--atac', type=click.Path(exists=True, resolve_path=True), required=True, help='ATAC.')
@click.option('-g', '--gtf', type=click.Path(exists=True, resolve_path=True), required=True, help='GTF.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('--rnabatchkey', type=click.STRING, default='sampleid', show_default=True, help='Sample batch key for RNA.')
@click.option('--rnahvgntop', type=click.INT, default=10000, show_default=True, help='Number of HVGs for RNA.')
@click.option('--rnahvgflavor', type=click.Choice(['seurat', 'cell_ranger', 'seurat_v3']), is_flag=False, flag_value='seurat_v3', default='seurat', show_default=True, help='HVG algorithm for RNA.')
@click.option('-c', '--atacpeakfrac', type=click.FLOAT, default=0.2, show_default=True, help='Fraction of variable peaks for ATAC.')
@click.option('--gtfby', type=click.STRING, default='gene_name', show_default=True, help='GTF field for gene name.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-R', '--refreshrna', is_flag=True, help='Process RNA.')
@click.option('-A', '--refreshatac', is_flag=True, help='Process ATAC.')
@click.option('-G', '--refreshgraph', is_flag=True, help='Process Graph guidance.')
def main(outdir, bname, rna, atac, gtf, seed, rnabatchkey, rnahvgntop, rnahvgflavor, atacpeakfrac, gtfby, condaenv, refreshrna, refreshatac, refreshgraph):
	"""
To preprocess RNA and ATAC for scglue co-embedding. RNA and ATAC .h5ad files are provided by "-r|--rna" and "-a|--atac".

\b
Example:
  rna=/project/rchen/shared_folder/HRCA/HRCA_BC.h5ad
  atac=/project/rchen/u236923/mrrdir/mwe/GLUEmwe/example/BC/preproc/ATAC/scrnascanpyconcat2h5ad/HRCA_BC_ATAC.h5ad
  gtf=/project/rchen/u236923/mrrdir/mwe/GLUEmwe/example/gtf/genes.gtf
  bname=HRCA_BC
  outdir=$(mrrdir.sh)
  slurmmab.sh -t 2 -m 200G -n gpu-8-1:0 -- gluernaatac2preproc -e scglue -d "$outdir" -b "$bname" -r "$rna" -a "$atac" -g "$gtf" --rnabatchkey sampleid --rnahvgntop 10000 --rnahvgflavor seurat -c 0.2 --gtfby gene_name

\b
  # An example to refresh RNA and Graph guidance but reuse ATAC
  rna=/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/ONONH/scglue/rnaraw_s10k/scrnah5adfiles2scviwkfl/ONONH_majorclass_rawCount.h5ad
  atac=/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/ONONH/scglue/gluernaatac2coembedwkfl/preproc_atac/ONONH_snATAC.h5ad
  gtf=/dfs3b/ruic20_lab/jinl14/resource/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz
  bname=ONONH
  outdir=$(mrrdir.sh)
  gluernaatac2preproc -e scglue -d "$outdir" -b "$bname" -r "$rna" -a "$atac" -g "$gtf" --rnabatchkey sampleid --rnahvgntop 10000 --rnahvgflavor seurat -c 0.2 --gtfby gene_name -R -G

\b
Note:
  1. RNA used "--rnabatchkey" for HVGs, while ATAC didn't consider sample batches for feature selection. (TODO)
  2. Refresh options (-R, -A, and -G) are used to redo the processing while reusing previously processed files for other modality.

\b
See also:
  Upstream:
    gluernaatac2preproc
  Depends:
    Python/scglue
    Python/Scanpy
    Python/epiScanpy

\b
Date: 2024/05/24
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/python/{scriptname}.py'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"rnafile='{rna}'",
		f"atacfile='{atac}'",
		f"gtffile='{gtf}'",
		f"seed={seed}",
		f"rnabatchkey='{rnabatchkey}'",
		f"rnahvgntop={rnahvgntop}",
		f"rnahvgflavor='{rnahvgflavor}'",
		f"atacpeakfrac={atacpeakfrac}",
		f"gtfby='{gtfby}'",
		f"refreshrna={refreshrna}",
		f"refreshatac={refreshatac}",
		f"refreshgraph={refreshgraph}",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
