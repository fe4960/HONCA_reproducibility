#!/bin/sh
source ~/.condainit
conda activate archr
R < HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/4heatmap_glia_age.R --no-save > HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/4heatmap_glia_age.out 2> HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/4heatmap_glia_age.err
