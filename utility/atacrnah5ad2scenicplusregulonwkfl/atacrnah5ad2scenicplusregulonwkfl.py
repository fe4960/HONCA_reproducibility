#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
import shutil
from exeBash.exeBash import exeBash
from pathlib import Path
import click
import datetime

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-c', '--configfile', type=click.Path(exists=False, resolve_path=True), help='Configuration file in the YAML format.')
@click.option('-w', '--workprofile', type=click.Path(exists=False, resolve_path=True), help='Snakemake workflow profile folder, e.g., for Slurm. E.g., ~/.config/snakemake/scenicplusslurm/.')
@click.option('-r', '--rna', type=click.Path(exists=False, resolve_path=True), required=True, help='RNA.')
@click.option('-t', '--numthreads', type=click.INT, default=12, show_default=True, help='Number of threads.')
@click.option('-n', '--dryrun', is_flag=True, help='Dry-run.')
@click.argument('atac', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, configfile, workprofile, rna, numthreads, dryrun, atac):
	"""
To calculate regulons from ATAC union peak and RNA reference.

ATAC: ATAC union peak .h5ad. See archrpeaksetfolder2unionpeakh5adwkfl.

\b
Example:
  atac=/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/ONONH/coembed/oligodendrocyte/postproc/peakset/archrproject2addpeakset/archrpeaksetfolder2unionpeakh5adwkfl/scrnah5adaddobsmfromtxt/ON_oligodendrocyte.h5ad
  rna=/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/ONONH/coembed/oligodendrocyte/rna/scrnah5ad2clean/ON_oligodendrocyte.h5ad
  bname=$(basename "$atac" .h5ad)
  outdir=$(mrrdir.sh)
  atacrnah5ad2scenicplusregulonwkfl -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -w ~/.config/snakemake/scenicplusslurm -t 12 -r "$rna" -- "$atac"

\b
Note:
  1. .h5ad file will be used for scGLUE.
  2. Parallel bug related to a the ray pacakge. ataccistopicobject2regionsets and ataccistopic2runscenicpluspipeline can not run inside a snakemake.
  2.1 Step 1. Run atacpeakset2cistopicmodel in parallel through slurm workprofile.

\b
  ```
  slurmucihpc.sh -t 2 -m 120G -w 3-00:00:00 -p free -- atacrnah5ad2scenicplusregulonwkfl -c config.yaml -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -t 40 -r "$rna" -w ~/.config/snakemake/scenicplusslurm -- "$atac"
  ```

\b
  2.2 Step 2. Run ataccistopicobject2regionsets and ataccistopic2runscenicpluspipeline on a single node without slurm workprofile.

\b
  ```
  slurmucihpc.sh -t 2 -m 420G -w 3-00:00:00 -p free -- atacrnah5ad2scenicplusregulonwkfl -c config.yaml -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -t 40 -r "$rna" -- "$atac"
  ```

\b
See also:
  Upstream:
    archrpeaksetfolder2unionpeakh5adwkfl
  Depends:
    Python/scenicplus

\b
Date: 2026/02/21
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	shutil.copytree(absdir / 'snakemake', outdir, symlinks=True, dirs_exist_ok=True, ignore=shutil.ignore_patterns('__pycache__'))
	if workprofile:
		workprofile=Path(workprofile).expanduser()
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')

	# common command string
	cmdstr=[
		f"snakemake",
		' '.join((
			"--config",
			f"rna='{rna}'",
			f"atac='{atac}'",
			f"outdir='{outdir}'",
			f"bname='{bname}'",
			f"condaenv='{condaenv}'",
			f"nowtimestr='{nowtimestr}'",
			)),
		f"--directory '{outdir}'",
		f"--jobs {numthreads}",
		f"--cores {numthreads}",
		*([f"--workflow-profile '{workprofile}'"] if workprofile and workprofile.exists() else []),
		]

	if configfile:
		cmdstr+=[f"--configfile '{configfile}'"]

	if dryrun: # dry-run
		cmdstr+=[f"--dry-run --printshellcmds"]

	else: # running
		cmdstr+=[
			f"--printshellcmds --debug-dag --skip-script-cleanup --software-deployment-method conda --keep-incomplete --keep-going",
			]

	os.chdir(outdir)
	return exeBash.callback(cmdstr=cmdstr, condaenv=None, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
