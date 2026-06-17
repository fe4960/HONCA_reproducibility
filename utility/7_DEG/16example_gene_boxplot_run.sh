#!/bin/sh

comd=HCA_ON/scripts/7_DEG/16example_gene_boxplot


mkdir $comd

c="Oligodendrocyte"

od="genexp_donor_subclass_final_update"

file=HCA_ON/data/7_DEG/${od}/ct_list_${c}

while read line
do
ct=$line

for gene in PIEZO2 YAZ1 TAFAZZIN ERG5 CCN2 EGR2 AKT1 MAPK1 PIK3CA PTEN MAP2K7 MAPK14 MAPK3  
do
#ct="OLIGO2_LRRC7+"
sbatch --account=ruic20_lab -p free --nodelist=hpc3-14-12  --mem=5GB --time=0-1 --error=${comd}/${line}_${gene}.err --output=${comd}/${line}_${gene}.out ${comd}.sh $c $ct $od $gene
done

done < $file
