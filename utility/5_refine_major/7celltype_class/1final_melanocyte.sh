#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution"

mkdir -p $err

#for cell in AC Astrocyte BC Fibroblast MG Macrophage Melanocyte Microglia Oligodendrocyte RGC RPE  
#for cell in Astrocyte
#for cell in Fibroblast
#for cell in Oligodendrocyte
for cell in Melanocyte
#for cell in Macrophage
do	
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean
indir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/

for res in  0.2
#for res in 0.2
do

h5ad=${indir}/${cell}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad

bname=${cell}_res_${res}_clean

#outdir=$indir
rmlist=rm_cluster_rs #"none" #rm_cluster_res_1
mk="sanes_mk"
save="t"
sbatch --mem=40GB -p standard --account=ruic20_lab --time=1-0 --output=${err}/${cell}_${res}.out --error=${err}/${cell}_${res}.err  HCA_ON/scripts/5_refine_major/5refine_resolution.sh $h5ad $res ${batch_key} $bname $outdir $rmlist $mk $save $indir
done
done
