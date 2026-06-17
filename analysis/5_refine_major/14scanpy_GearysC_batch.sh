#!/bin/sh
err=HCA_ON/scripts/5_refine_major/14scanpy_GearysC
mkdir $err
comd=HCA_ON/scripts/5_refine_major/14scanpy_GearysC.sh
#for cell in  Astrocyte Cone HC MG Macrophage Mural_cell Schwann_cell  NK_T_cell  Oligodendrocyte_precursor_cell Microglia Endothelial_cell AC BC RGC RPE Fibroblast  Oligodendrocyte
#for cell in Oligodendrocyte
#for cell in Melanocyte
for cell in major
do
sbatch --mem=200GB -p standard --account=ruic20_lab --time=0-2 --output=${err}/${cell}.out --error=${err}/${cell}.err  $comd $cell
done
