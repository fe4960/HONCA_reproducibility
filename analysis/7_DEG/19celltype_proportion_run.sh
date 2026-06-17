#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/19celltype_proportion
mkdir $comd
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor

for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell
do
sh sb.sh -c ${comd}.sh -e $comd -j $c -- $c $meta
done
