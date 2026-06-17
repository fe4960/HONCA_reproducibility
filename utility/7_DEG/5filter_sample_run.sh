#!/bin/sh

comd=HCA_ON/scripts/7_DEG/5filter_sample

mkdir $comd
od="genexp_donor_subclass_final_new"

#od="genexp_donor_subclass_final_update"
#for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell #  Endothelial_cell  Oligodendrocyte Fibroblast 
#for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell Endothelial_cell Fibroblast Microglia RGC MG  Macrophage
for c in Astrocyte #Oligodendrocyte Fibroblast #  Endothelial_cell  Oligodendrocyte Fibroblast 
#for c in Endothelial_cell
#for c in major
do

#file=HCA_ON/data/7_DEG/ct_list_${c}
#dir=HCA_ON/data/7_DEG/${c}_majorclass/
#dir=HCA_ON/data/7_DEG/${c}_majorclass_sample/

file=HCA_ON/data/7_DEG/${od}/ct_list_${c}2
dir=HCA_ON/data/7_DEG/${od}/${c}_subclass/

	
#file=HCA_ON/data/7_DEG/genexp_donor_subclass/ct_list_${c}
#dir=HCA_ON/data/7_DEG/genexp_donor_subclass/${c}_subclass/


#file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${c}
#dir=HCA_ON/data/7_DEG/genexp_sample_subclass/${c}_subclass/

#while read line
#do
#ct=$line
sbatch -p free --mem=6GB --time=0-1 --error=${comd}/${c}.err --output=${comd}/${c}.out ${comd}.sh $file $dir
#done < $file
done
