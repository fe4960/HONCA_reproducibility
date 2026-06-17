#!/bin/sh
source ~/.condainit
conda activate scvi

hvg=2000
max_epoch="none"
fl="seurat_v3"
res=0.3
main="/dfs3b/ruic20_lab/junw42/"
file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
err="/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/6refine_scvi_plot"
mkdir $err
#for celltype in Endothelial_cell
for celltype in Fibroblast
#for celltype in Astrocyte
#while read celltype
#for celltype in BC  Melanocyte Oligodendrocyte_precursor_cell
#for celltype in Macrophage  #AC RPE RGC Oligodendrocyte_precursor_cell MG Endothelial_cell BC
do
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean"
bname="${celltype}_hvg10000_epochnone_seurat_v3_rs_1_span_1_subclass"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mk="sanes_mk"
sbatch -p standard --mem=30GB --time=0-2 --error=${err}/${bname}.err --output=${err}/${bname}.out --account=ruic20_lab /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/6refine_scvi_plot.sh $indir $mk $outdir $bname
done
#done < $file
