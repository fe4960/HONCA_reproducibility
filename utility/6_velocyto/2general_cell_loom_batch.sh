#!/bin/sh

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/2general_cell_loom
mkdir $err
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/2general_cell_loom.sh

#indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
#h5ad="Oligodendrocyte_subclass"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir

#indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
#h5ad="Oligodendrocyte_precursor_cell_subclass"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir



indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
h5ad="Macrophage_subclass_sb"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir



indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/MG/clean/"
h5ad="MG_subclass_sb"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
#h5ad="Microglia_subclass_sb"
h5ad="Microglia_subclass"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
h5ad="Astrocyte_subclass_new"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/"
h5ad="RPE_subclass"

#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Schwann_cell/clean/"
h5ad="Schwann_cell_subclass"


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
#h5ad="Astrocyte_hvg_2000_fl_seurat_v3_res_0.4_sd_7_ononh_retina_nonsb_subclass"
h5ad="Astrocyte_subclass_new_update"
#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
h5ad="Oligodendrocyte_subclass_seurat"

#sbatch -p free-gpu --mem=100GB --time=1-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $h5ad $indir

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
h5ad="Oligodendrocyte_precursor_cell_subclass_seurat"

#sbatch -p free-gpu --mem=100GB --time=1-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $h5ad $indir


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
h5ad="oligo_opc_seurat_new_scvi"

sbatch -p free-gpu --mem=100GB --time=1-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $h5ad $indir

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
h5ad="Astrocyte_subclass_new5_clean"

#sbatch -p free-gpu --mem=100GB --time=0-20 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd  $h5ad $indir


