#!/bin/sh

comd=HCA_ON/scripts/7_DEG/6DEG_dream_donor_sample

mkdir $comd
#for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell Endothelial_cell Fibroblast Microglia RGC MG  Macrophage      	
#do
#od="genexp_donor_subclass_final_update"

c="major"

file=HCA_ON/data/7_DEG/ct_list_major

for line in Astrocyte Microglia Fibroblast Macrophage  
#while read line
do
ct=$line
#sbatch --account=ruic20_lab -p standard  --nodelist=hpc3-14-06   --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od
#sbatch --account=ruic20_lab -p standard  --nodelist=hpc3-21-15   --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od
sbatch --account=ruic20_lab -p free-gpu     --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

#done < $file
done



