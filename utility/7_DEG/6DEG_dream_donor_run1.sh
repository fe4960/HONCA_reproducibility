#!/bin/sh

comd=HCA_ON/scripts/7_DEG/6DEG_dream_donor

mkdir $comd
#for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell Endothelial_cell Fibroblast Microglia RGC MG  Macrophage      	
#do
od="genexp_donor_subclass_final_update"

c="Astrocyte"

file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_0

while read line
do
ct=$line
sbatch --account=ruic20_lab -p standard  --nodelist=hpc3-14-12   --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file
#done


c="Macrophage"

file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_12

while read line
do
ct=$line
#sbatch --account=ruic20_lab -p standard --nodelist=hpc3-14-03 --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file

c="Oligodendrocyte_precursor_cell"

file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_13

while read line
do
ct=$line
######sbatch --account=ruic20_lab -p standard --nodelist=hpc3-14-12 --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file

c="Microglia"

file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_14

while read line
do
ct=$line
###########sbatch --account=ruic20_lab -p standard --nodelist=hpc3-20-26 --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file

c="Oligodendrocyte"

file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_15

while read line
do
ct=$line
sbatch --account=ruic20_lab -p standard --nodelist=hpc3-14-12 --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file

c="RGC"

file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_16

while read line
do
ct=$line
########sbatch --account=ruic20_lab -p standard --nodelist=hpc3-19-17 --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file

c="Endothelial_cell"


file=HCA_ON/data/7_DEG/${od}/ct_list_all_4_17

while read line
do
ct=$line
#sbatch --account=ruic20_lab -p standard --nodelist=hpc3-14-12   --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

done < $file


