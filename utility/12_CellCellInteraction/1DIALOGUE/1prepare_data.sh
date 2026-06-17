#!/bin/sh
source ~/.condainit
conda activate scvi
dir="/dfs3b/ruic20_lab/junw42"
#fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_final"
#fn=${dir}/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_rawcount
fn=${dir}/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_Astrocyte_rawcount
label="celltype"
#label="majorclass"
cellnum=20000 #5000
outdir=${dir}/HCA_ON/data/10_DIALOGUE/
mkdir -p $outdir
#ct="ONONH_majorclass_snRNA"
#rml=${outdir}/rm_list
rml="none"
#out=${outdir}/ONONH_${label}_${cellnum}
out=${outdir}/ONONH_Astrocyte_${label}_${cellnum}

ct="n"
python /dfs3b/ruic20_lab/junw42/human_ret_anc/scripts/snRNA/7downsample.py  $fn $label $cellnum $out $rml   $ct
