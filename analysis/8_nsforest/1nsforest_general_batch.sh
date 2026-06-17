#!/bin/sh

comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/8_nsforest/1nsforest_general.sh
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/8_nsforest/1nsforest_general
mkdir $err

#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/"
#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
####fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
#fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"

#fn="ONONH_all_normcount"
#######fn="Fibroblast_subclass_LC_other1"
#fn=Oligodendrocyte_subclass_new
#fn="Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg_update_full"
#fn="Oligodendrocyte_precursor_cell_subclass_seurat"
#fn="Astrocyte_subclass_new5_clean"

#fn="Astrocyte_subclass_new4_clean"
#fn="Fibroblast_subclass_clean_new"
#fn=Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_full_clean_update
#fn="Oligodendrocyte_subclass_seurat"
#output="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/nsforest"

mkdir -p $output


#############cluster="majorclass"

#########sbatch --account=ruic20_lab -p standard --time=0-12 --mem=200GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 


########cluster="celltype"

#for 
fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/"
fn="Endothelial_cell_subclass_sb_seurat_rmRPE"
cluster="subclass"
output=${fold}/nsforest

##########sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=20GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 



fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/"
fn="Mural_cell_subclass_rmRPE"
cluster="subclass"
output=${fold}/nsforest

###########sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=20GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 


fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
fn="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location1"
#fn="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_recluster"
cluster="subclass"
output=${fold}/nsforest

#cluster="leiden1"
sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=50GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 


fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
fn="Oligodendrocyte_precursor_cell_subclass_seurat_cycling_rmRPE"
cluster="subclass"
output=${fold}/nsforest

###############sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=40GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 



fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
fn="Oligodendrocyte_subclass_seurat"
cluster="subclass"
output=${fold}/nsforest

############3sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=70GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 



fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
fn="Macrophage_subclass_sb_clean_rmRPE"
cluster="subclass"
output=${fold}/nsforest

############sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=20GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 



fold="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
fn="Microglia_subclass_sb_clean"
cluster="subclass"
output=${fold}/nsforest

############sbatch --account=ruic20_lab  -p standard --time=0-4 --mem=20GB --error=${err}/${cluster}.err --output=${err}/${cluster}.out $comd $fold $fn $output $cluster 

