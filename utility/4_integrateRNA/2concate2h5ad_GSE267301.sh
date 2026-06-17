#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"
#bname="GSE167494"
bname="GSE167494_more"

outdir=${main}/HCA_ON/data/9_Brain_h5ad/${bname}

#outdir=${main}/HCA_ON/data/scvi/RNA/raw/
#bname="GSE267301"
mkdir -p $outdir
#samplelist=${main}/HCA_ON/data/sample_list/RNA_sample_list_meta.gz
#samplelist=${main}/HCA_ON/data/2_cellqc/pub/GSE167494_h5ad_cellqc
samplelist=${main}/HCA_ON/data/2_cellqc/pub/GSE167494_h5ad_cellqc2_sort

#samplelist=${main}/HCA_ON/data/2_cellqc/pub/GSE267301_h5ad
#samplelist=${main}/HCA_ON/data/2_cellqc/pub/GSE267301_h5ad

python ${main}/HCA_ON/scripts/4_integrateRNA/2concate2h5ad.py $samplelist $bname $outdir


#scrnascanpyconcat2h5ad.sh -t 12 -d "$outdir" -b "$bname" -- $fname
