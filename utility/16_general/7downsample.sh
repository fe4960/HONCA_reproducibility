#!/bin/sh
source ~/.condainit
conda activate scvi

#fn=sys.argv[1]
#cluster_key=sys.argv[2]
#target_cells = int(sys.argv[3])
#out=sys.argv[4]

#label=sys.argv[5]
#celltype=sys.argv[6]


python HCA_ON/scripts/16_general/7downsample.py $1 $2 $3 $4 $5 $6
