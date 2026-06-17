#!/bin/sh
comd=HCA_ON/scripts/7_DEG/12GO_analysis_mash_ct
mkdir $comd
#dem="age"
dem="age_cpm1"
#dem="age_range_30_60_cpm5"

#for cell in Astrocyte Fibroblast Oligodendrocyte Endothelial_cell
for cell in major
do
file=HCA_ON/data/7_DEG/ct_list_${cell}

#file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${cell}
j=1
while read line
do
#j=$(($j+1))
ct=$line
sbatch --nodelist=hpc3-l18-01 --account=ruic20_lab -p standard --mem=10GB --time=0-3 --error=${comd}/${cell}_${ct}.err --output=${comd}/${cell}_${ct}.out ${comd}.sh $cell $ct $j $dem
done < $file
done
