#!/bin/sh
err=HCA_ON/scripts/5_refine_major/13subclass_umap
comd=${err}.sh
mkdir $err

#for cell in Oligodendrocyte Fibroblast Microglia Endothelial_cell BC Oligodendrocyte_precursor_cell Mural_cell MG AC RPE Macrophage HC Cone NK_T_cell RGC Schwann_cell  Melanocyte B_cell 
for cell in Endothelial_cell
do
sbatch  --account=ruic20_lab -p free --mem=50GB --time=0-1 --error=${err}/${cell}.err --output=${err}/${cell}.out $comd $cell
done
