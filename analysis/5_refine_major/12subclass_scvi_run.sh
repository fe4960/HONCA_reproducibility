#!/bin/sh
err=HCA_ON/scripts/5_refine_major/12subclass_scvi
comd=${err}.sh
mkdir $err

#for cell in Oligodendrocyte Fibroblast Microglia Endothelial_cell BC Oligodendrocyte_precursor_cell Mural_cell MG AC RPE Macrophage HC Cone NK_T_cell RGC Schwann_cell  Melanocyte B_cell  
#for cell in Astrocyte
for cell in Oligodendrocyte
do
sbatch --gres=gpu:1 --account=ruic20_lab_gpu -p free-gpu --mem=50GB --time=0-2 --error=${err}/${cell}.err --output=${err}/${cell}.out $comd $cell
done
