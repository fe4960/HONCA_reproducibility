#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution_dominant_sample"

mkdir -p $err

for cell in AC Astrocyte BC Fibroblast MG Macrophage Melanocyte Microglia Oligodendrocyte RGC RPE  
do	
indir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/

for res in 1 #0.1 0.2 0.3
do

#h5ad=${indir}/${cell}_hvg2k_epoch20_scvi_trg.h5ad

#bname=${cell}_res_${res}
bname=${cell}_hvg2k_epoch20_${cell}_scvi
outdir=$indir

sbatch --mem=40GB -p standard --account=ruic20_lab --time=4-0 --output=${err}/${cell}_${res}.out --error=${err}/${cell}_${res}.err  HCA_ON/scripts/5_refine_major/5refine_resolution_dominant_sample.sh  $bname $outdir
done
done
