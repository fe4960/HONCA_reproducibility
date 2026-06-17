#!/bin/sh
source ~/.condainit
conda activate nsforest
#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/"
#fn="Endothelial_cell_subclass_sb"
#####fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
######fn="Fibroblast_subclass_LC_other"

fold=$1

fn=$2

output=$3

cluster=$4 

#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/"
#fn="ONONH_all_normcount"
#output="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/nsforest"
mkdir -p $output
#cluster="majorclass"

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/8_nsforest/1nsforest_general.py $fold $fn $output $cluster
