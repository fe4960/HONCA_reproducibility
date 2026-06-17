#!/bin/sh
source ~/.condainit
conda activate qtl
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor

mkdir $comd

main_cell=$1
cell=$2
od=$3
dirInName=${cell}_dream
fn=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor
rml=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${od}/${main_cell}_subclass/exp_${cell}_cor_Norm_all_flt
dirIn=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${od}/${main_cell}_subclass_DEG/${cell}/
indir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/${od}/${main_cell}_subclass


#echo "Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor.R $cell $dirInName $fn $rml $dirIn $indir "

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor.R $cell $dirInName $fn $rml $dirIn $indir

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/6DEG_dream_donor_age_range.R $cell $dirInName $fn $rml $dirIn $indir

R --vanilla --args $cell $dirInName $fn $rml $dirIn $indir < ${comd}.R > ${comd}/${cell}.out 2> ${comd}/${cell}.err 

