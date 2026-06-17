#!/bin/sh

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/3scvelo_general_batch
mkdir $err
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/3scvelo_general.sh

#hvg=5000
hvg=2000

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
h5ad="Oligodendrocyte_subclass"

#sbatch -p standard --mem=200GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw $hvg

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
h5ad="Oligodendrocyte_precursor_cell_subclass"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad  ${h5ad}_hvg_${hvg}_raw $hvg 



indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
h5ad="Macrophage_subclass_sb"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/MG/clean/"
h5ad="MG_subclass_sb"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg



indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
h5ad="Microglia_subclass_sb"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
#h5ad="Microglia_subclass_sb"
h5ad="Microglia_subclass"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
h5ad="Astrocyte_subclass"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/"
h5ad="RPE_subclass"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Schwann_cell/clean/"
#h5ad="Schwann_cell_subclass_sb"
h5ad="Schwann_cell_subclass"

#sbatch -p standard --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
#h5ad="Schwann_cell_subclass_sb"
h5ad="Astrocyte_subclass_new_update"

sbatch -p standard --mem=100GB --time=4-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $indir $h5ad ${h5ad}_hvg_${hvg}_raw  $hvg
