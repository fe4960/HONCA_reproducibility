#!/bin/sh
comd=HCA_ON/scripts/20_hotspot/4GO_enrichR
mkdir $comd
#cell="Oligodendrocyte_precursor_cell_permut"
#cell="Oligodendrocyte_precursor_cell_fullcell_10kgene"
#cell="Astrocyte"
#cell="Oligodendrocyte_fullcell_5kgene"
bg1="t"
#for ct in 100_hvg_5000 160_hvg_5000 200_hvg_5000 250_hvg_5000
for cell in  Astrocyte  # OLIGO2_LRRC7+_permute_metacell_200 OLIGO1_permute_metacell_200 OLIGO1_SVEP1+_permute_metacell_200 OLIGO1_RBFOX1+_permute_metacell_200
do
#t="age"
t="scvi"
for ct in  100_hvg_5000 #   350_hvg_5000   300_hvg_5000   250_hvg_5000 200_hvg_5000   160_hvg_5000  120_hvg_5000   100_hvg_5000 #  160_hvg_10000 250_hvg_10000 350_hvg_10000 
do
sbatch  --account=ruic20_lab -p free --mem=10GB --time=0-1 --error=${comd}/${cell}_${ct}.err --output=${comd}/${cell}_${ct}.out ${comd}.sh $cell $ct $bg1 $t
done
done

