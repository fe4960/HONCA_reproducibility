#!/bin/sh
comd=HCA_ON/scripts/10_compare_SEAD/6data_source_general_new

for c in Astrocyte Microglia Oligodendrocyte Oligodendrocyte_precursor_cell
do
sbatch -p free --mem=20GB --time=0-1 --output=${comd}/${c}.out --error=${comd}/${c}.err ${comd}.sh $c
done
