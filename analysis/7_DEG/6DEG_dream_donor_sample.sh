#!/bin/sh
source ~/.condainit
conda activate qtl

main_cell=$1
cell=$2
dirInName=${cell}_dream_major_sample
fn=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid
rml=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${main_cell}_majorclass_sample/exp_${cell}_cor_Norm_all_flt
dirIn=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${main_cell}_majorclass_sample_DEG/${cell}/
indir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${main_cell}_majorclass_sample/


#echo "Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream.R $cell $dirInName $fn $rml $dirIn $indir "

Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor_sample.R $cell $dirInName $fn $rml $dirIn $indir

