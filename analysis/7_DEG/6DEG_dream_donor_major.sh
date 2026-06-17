#!/bin/sh
source ~/.condainit
conda activate qtl

main_cell=$1
cell=$2
od=$3
dirInName=${cell}_dream
fn=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor
rml=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${main_cell}_majorclass/exp_${cell}_cor_Norm_all_flt
dirIn=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${main_cell}_majorclass_DEG/${cell}/
indir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${main_cell}_majorclass

mkdir -p $dirIn

#echo "Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor.R $cell $dirInName $fn $rml $dirIn $indir "

Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor.R $cell $dirInName $fn $rml $dirIn $indir

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor_age_range.R $cell $dirInName $fn $rml $dirIn $indir

