#!/bin/sh
  
source ~/.condainit
conda activate scvi

cell=$1

python HCA_ON/scripts/5_refine_major/14scanpy_GearysC.py $cell
