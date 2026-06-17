#!/bin/sh
source ~/.condainit
conda activate veloc


python HCA_ON/scripts/6_velocyto/3scvelo_general_plot_astro.py 

#$indir $h5ad $bname $hvg 

#python HCA_ON/scripts/6_velocyto/3scvelo_general_oligo_opc.py $indir $h5ad $bname $hvg 
