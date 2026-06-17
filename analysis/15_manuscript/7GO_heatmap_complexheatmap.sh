#!/bin/sh

source ~/.condainit 
conda activate archr
err=HCA_ON/scripts/15_manuscript/7GO_heatmap_complexheatmap
mkdir $err
for c in Oligodendrocyte Oligodendrocyte_precursor_cell Fibroblast Astrocyte
do
fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${c}/clean/${c}_subclass_seurat_trg500_enrichR/${c}_GO_top15_MSigDB_Hallmark_2020_up"

#fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${c}/clean/${c}_subclass_seurat_trg500_enrichR/${c}_GO_top10_GO_Biological_Process_2023_up"

#fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${c}/clean/${c}_subclass_seurat_trg500_enrichR/${c}_GO_top10_MSigDB_Hallmark_2020_up"
Rscript --vanilla ${err}.R $fn > ${err}/${c}.out 2> ${err}/${c}.err
done
