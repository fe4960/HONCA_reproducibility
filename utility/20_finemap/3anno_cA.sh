#!/bin/sh

source ~/.condainit

conda activate archr

R < HCA_ON/scripts/20_finemap/3anno_cA.R  --no-save > HCA_ON/scripts/20_finemap/3anno_cA.out 2> HCA_ON/scripts/20_finemap/3anno_cA.err
