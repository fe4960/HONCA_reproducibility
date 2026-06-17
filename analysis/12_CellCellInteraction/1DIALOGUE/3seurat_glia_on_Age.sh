#!/bin/sh
source ~/.condainit
conda activate archr
comd="/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/3seurat_glia_on_Age"
R < ${comd}.R --no-save > ${comd}.out 2> ${comd}.err 
