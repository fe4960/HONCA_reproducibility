#!/bin/sh

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/4scvelo_plot_batch
mkdir $err
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/4scvelo_plot.sh

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
h5ad="Oligodendrocyte_subclass"

sbatch -p standard --mem=200GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir

h5ad="oligo_opc"

sbatch -p standard --mem=200GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
h5ad="Oligodendrocyte_precursor_cell_subclass"

sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir



indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
h5ad="Macrophage_subclass"

sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir

