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
@click.option('-g', '--graphguidance', type=click.Path(exists=True, resolve_path=True), required=True, help='Graph guidance.')
@click.option('-s', '--seed', type=click.INT, default=12345, show_default=True, help='Random seed.')
@click.option('-l', '--rnalabel', type=click.STRING, default='celltype', show_default=True, help='Cell label for RNA.')
@click.option('-L', '--ataclabel', type=click.STRING, default='celltype', show_default=True, help='Cell label for ATAC.')
@click.option('-R', '--rnabatchkey', type=click.STRING, default='sampleid', show_default=True, help='Sample batch key for RNA.')
@click.option('-A', '--atacbatchkey', type=click.STRING, default='sampleid', show_default=True, help='Sample batch key for ATAC.')
@click.option('-m', '--rnareps', type=click.STRING, multiple=True, default=['X_scVI', 'X_pca'], show_default=True, help='RNA representation.')
@click.option('-n', '--atacreps', type=click.STRING, multiple=True, default=['Harmony', 'IterativeLSI', 'X_lsi'], show_default=True, help='ATAC representation.')
@click.option('--latent_dim', type=click.INT, default=30, show_default=True, help='Latent dimensionality.')
@click.option('--h_depth', type=click.INT, default=2, show_default=True, help='Hidden layer depth for encoder and discriminator.')
@click.option('--data_batch_size', type=click.INT, default=128, show_default=True, help='Number of cells in each data minibatch.')
@click.option('--graph_batch_size', type=click.INT, default=128, show_default=True, help='Number of edges in each graph minibatch.')
@click.option('--align_burnin', type=click.INT, default=32, show_default=True, help='Number of epochs to wait before starting alignment.')
@click.option('--max_epochs', type=click.INT, default=200, show_default=True, help='Maximal number of epochs.')
@click.option('--patience', type=click.INT, default=16, show_default=True, help='Patience of early stopping.')
@click.option('--reduce_lr_patience', type=click.INT, default=8, show_default=True, help='Patience to reduce learning rate.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
def main(
		outdir,
		bname,
		rna,
		atac,
		graphguidance,
		seed,
		rnalabel,
		ataclabel,
		rnabatchkey,
		atacbatchkey,
		rnareps,
		atacreps,
		latent_dim,
		h_depth,
		data_batch_size,
		graph_batch_size,
		align_burnin,
		max_epochs,
		patience,
		reduce_lr_patience,
		numthreads,
		condaenv
		):
	"""
To construct scglue co-embedding from preprocessed files. See gluernaatac2preproc for preprocessing.

\b
Example:
  indir=/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/ONONH/scglue/gluernaatac2coembedwkfl
  rna=$indir/preproc_rna/ONONH_snRNA.h5ad
  atac=$indir/preproc_atac/ONONH_snATAC.h5ad
  graphguidance=$indir/preproc_graphguidance/ONONH_snRNA_ONONH_snATAC.graphml.gz
  bname=ONONH_glue
  outdir=$(mrrdir.sh)
  function cmd {
  local rnarep=$1
  local atacrep=$2
  slurmucihpc.sh -t 2 -m 380G -w 48:00:00 -p free-gpu -e jinl14 -- glueh5adgraphatacrnatype2coembed -e scglue -d "$outdir" -b "$bname" -r "$rna" -a "$atac" -g "$graphguidance" -l majorclass -L majorclass -R sampleid -A sampleid -m "$rnarep" -n "$atacrep"
  }
  
  source env_parallel.bash
  env_parallel cmd ::: X_scVI X_pca ::: Harmony IterativeLSI X_lsi

\b
Note:
  1. This wrapper is a helper wrapper for generating co-embedding visualizations after annotating ATAC cells. To annotate ATAC cells, use glueh5adgraph2coembed.
  2. glueh5adgraph2coembed performed co-embedding in an unsupervised mode, while this wrapper uses supervised cell types.

\b
See also:
  Related:
    glueh5adgraph2coembed
  Upstream:
    gluernaatac2preproc
    gluernaatac2coembedwkfl
  Depends:
    Python/scglue

\b
Date: 2025/04/30
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
		f"guidefile='{graphguidance}'",
		f"seed={seed}",
		f"rnalabel='{rnalabel}'",
		f"ataclabel='{ataclabel}'",
		f"rnabatchkey='{rnabatchkey}'",
		f"atacbatchkey='{atacbatchkey}'",
		f"rnareps={rnareps}",
		f"atacreps={atacreps}",
		f"latent_dim={latent_dim}",
		f"h_depth={h_depth}",
		f"data_batch_size={data_batch_size}",
		f"graph_batch_size={graph_batch_size}",
		f"align_burnin={align_burnin}",
		f"max_epochs={max_epochs}",
		f"patience={patience}",
		f"reduce_lr_patience={reduce_lr_patience}",
		f"numthreads={numthreads}",
		]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)
	return exePython.callback(exprs, script=script, condaenv=condaenv, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
