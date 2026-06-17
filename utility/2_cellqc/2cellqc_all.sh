#!/bin/sh

source ~/.condainit
conda activate cellqc


#tmp="/pub/junw42/proj/HCA_ON/scripts/1cellqc/ON_chen_crpath_cellqc.txt_reform2"

#tmp="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/cellranger/cellranger_list_full"

#tmp="/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE267301/cellranger/cellranger_list_full"
#tmp="/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE167494/cellranger_sam/cellranger_list_full"
tmp="/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE167494/cellranger_sam/cellranger_list_full1"

#out="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellqc/pub/"
#out="/dfs3b/ruic20_lab/junw42/HCA_ON/data/2_cellqc/pub/GSE267301"
out="/dfs3b/ruic20_lab/junw42/HCA_ON/data/2_cellqc/pub/GSE167494"


mkdir -p $out
cellqc -d $out -t 10 -c /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellqc/ctrl.yaml  -- $tmp


