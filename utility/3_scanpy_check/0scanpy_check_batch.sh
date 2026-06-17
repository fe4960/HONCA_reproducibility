#!/bin/sh
source ~/.condainit
conda activate scglue
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/3_scanpy_check/0scanpy_check_batch.py $1 $2 $3 $4 

