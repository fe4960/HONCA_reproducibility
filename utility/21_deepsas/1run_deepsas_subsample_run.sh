#!/bin/sh

err=HCA_ON/scripts/21_deepsas/1run_deepsas_subsample


#ct="Oligodendrocyte"
#sh sb.sh -c ${err}.sh -p g -g t -m 60 -t 0-5 -a ruic20_lab_gpu -e $err -j $ct -- $ct



ct="Oligodendrocyte_precursor_cell"
#sh sb.sh -c ${err}.sh -p g -g t -m 60 -t 0-5 -a ruic20_lab_gpu -e $err -j $ct -- $ct


ct="Fibroblast"

#sh sb.sh -c ${err}.sh -p g -g t -m 70 -t 0-20 -a ruic20_lab_gpu -e $err -j $ct -- $ct


#ct="Astrocyte"

#for ct in Microglia Macrophage Mural_cell Endothelial_cell
for ct in Oligodendrocyte #Astrocyte
do
sh sb.sh -c ${err}.sh -p g -g t -m 70 -t 1-2 -a ruic20_lab_gpu -e $err -j $ct -- $ct

done
