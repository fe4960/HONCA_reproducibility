#!/bin/sh
source ~/.condainit
conda activate archr
R < HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/4heatmap_glia_age_on.R --no-save > HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/4heatmap_glia_age_on.out 2> HCA_ON/scripts/12_CellCellInteraction/1DIALOGUE/4heatmap_glia_age_on.err
