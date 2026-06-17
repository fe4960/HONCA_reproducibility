#!/bin/sh
comd=HCA_ON/scripts/7_DEG/14GO_enrichR_mash_ct
mkdir $comd
dem="age"
for cell in all #Astrocyte Fibroblast Oligodendrocyte Endothelial_cell
do
file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${cell}
j=0
while read line
do
j=$(($j+1))
ct=$line
sbatch  --account=ruic20_lab -p standard --mem=10GB --time=0-3 --error=${comd}/${cell}_${ct}.err --output=${comd}/${cell}_${ct}.out ${comd}.sh $cell $ct $j $dem

#sbatch --nodelist=hpc3-l18-01 --account=ruic20_lab -p standard --mem=10GB --time=0-3 --error=${comd}/${cell}_${ct}.err --output=${comd}/${cell}_${ct}.out ${comd}.sh $cell $ct $j $dem
done < $file
done
