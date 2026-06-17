#!/bin/sh

comd=HCA_ON/scripts/7_DEG/6DEG_dream_donor

mkdir $comd
###for c in Oligodendrocyte ###Oligodendrocyte_precursor_cell Astrocyte
####do
line="opc"
od="genexp_donor_subclass" #"genexp_donor_subclass_final"
file=HCA_ON/data/7_DEG/${od}/ct_list_${c}_1
c="all"
#####while read line
####do
ct=Oligodendrocyte_precursor_cell ##$line
sbatch --account=ruic20_lab -p free-gpu --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct $od

#####done < $file
#####done

