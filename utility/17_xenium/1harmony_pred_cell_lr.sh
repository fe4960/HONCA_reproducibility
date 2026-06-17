#!/bin/sh

source ~/.condainit
conda activate scvi

python HCA_ON/scripts/17_xenium/1harmony_pred_cell_lr.py
