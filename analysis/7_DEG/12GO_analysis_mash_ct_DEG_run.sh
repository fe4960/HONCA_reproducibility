#!/bin/sh
comd=HCA_ON/scripts/7_DEG/12GO_analysis_mash_ct_DEG
mkdir $comd
dem="age"
#dem="age_range_30_60_cpm5"
j=1
for dem in age_cpm1 age_range_30_60_cpm5
do
#for cell in Astrocyte Fibroblast Oligodendrocyte Endothelial_cell
for cell in major
do
file=HCA_ON/data/7_DEG/ct_list_${cell}

#file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${cell}
#j=1
while read line
do
#j=$(($j+1))
ct=$line

sbatch --nodelist=hpc3-l18-04 --account=ruic20_lab -p standard --mem=5GB --time=0-2 --error=${comd}/${cell}_${ct}_${dem}.err --output=${comd}/${cell}_${ct}_${dem}.out ${comd}.sh $cell $ct $j $dem


#sbatch --nodelist=hpc3-l18-01 --account=ruic20_lab -p standard --mem=10GB --time=0-3 --error=${comd}/${cell}_${ct}_${dem}.err --output=${comd}/${cell}_${ct}_${dem}.out ${comd}.sh $cell $ct $j $dem
done < $file
done

done
