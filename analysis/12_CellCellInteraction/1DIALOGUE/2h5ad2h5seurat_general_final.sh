#!/bin/sh
comd=HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/2h5ad2h5seurat_general_final.R
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/10_DIALOGUE/
#sample="ONONH_majorclass_5000"
#sample="ONONH_Astrocyte_celltype_5000"

source ~/.condainit
conda activate archr

#for sample in ONONH_Astrocyte_celltype_40000  ONONH_majorclass_40000
for sample in ONONH_Astrocyte_celltype_20000
#for sample in ONONH_Astrocyte_celltype_100000  ONONH_majorclass_100000
do
Rscript --vanilla $comd $sample $dir > HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/2h5ad2h5seurat_general_final/${sample}.out 2> HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/2h5ad2h5seurat_general_final/${sample}.err
done
