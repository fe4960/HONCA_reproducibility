#!/bin/sh

comd=HCA_ON/scripts/21_deepsas/2general_3_tables

mkdir $comd

ct="Microglia"
end=1
#sh sb.sh -c ${comd}.sh -p g -g t -t 0-4 -m 50 -a ruic20_lab_gpu -e $comd -j $ct --  $ct $end

ct="Endothelial_cell"

#sh sb.sh -c ${comd}.sh -p g -g t -t 0-4 -m 50 -a ruic20_lab_gpu -e $comd -j $ct --  $ct $end

end=3

ct="Fibroblast"

#sh sb.sh -c ${comd}.sh -p g -g t -t 0-4 -m 50 -a ruic20_lab_gpu -e $comd -j $ct --  $ct $end

ct="Mural_cell"

end=0
sh sb.sh -c ${comd}.sh -p g -g t -t 0-4 -m 50 -a ruic20_lab_gpu -e $comd -j $ct --  $ct $end

