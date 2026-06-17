#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"
outdir=${main}/HCA_ON/data/scvi/RNA/raw/
bname="ONONH_RNA_full"
mkdir -p $outdir
samplelist=${main}/HCA_ON/data/sample_list/RNA_sample_list_meta.gz
#samplelist=${main}/HCA_ON/data/2_cellqc/pub/GSE267301_h5ad

python ${main}/HCA_ON/scripts/4_integrateRNA/2concate2h5ad.py $samplelist $bname $outdir


#scrnascanpyconcat2h5ad.sh -t 12 -d "$outdir" -b "$bname" -- $fname
