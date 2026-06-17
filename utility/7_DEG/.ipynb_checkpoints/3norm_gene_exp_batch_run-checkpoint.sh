#!/bin/sh

comd=HCA_ON/scripts/7_DEG/3norm_gene_exp_batch

mkdir $comd

#for c in all Astrocyte  Endothelial_cell  Oligodendrocyte Fibroblast 
#for c in all Endothelial_cell Oligodendrocyte
#for c in Astrocyte #all #Endothelial_cell
#for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell
for c in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell Endothelial_cell Fibroblast Microglia RGC MG  Macrophage
do
#file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${c}
#dir=HCA_ON/data/7_DEG/genexp_sample_subclass/${c}_subclass/
file=HCA_ON/data/7_DEG/genexp_donor_subclass_final_update/ct_list_${c}
dir=HCA_ON/data/7_DEG/genexp_donor_subclass_final_update/${c}_subclass/

#while read line
#do
#ct=$line
sbatch -p free-gpu --mem=10GB --time=0-1 --error=${comd}/${c}.err --output=${comd}/${c}.out ${comd}.sh $file $dir
#done < $file
done
