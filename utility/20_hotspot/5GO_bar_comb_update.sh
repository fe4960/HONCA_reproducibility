#!/bin/sh

#source ~/.condainit

#conda activate qtl


source ~/.condainit
conda activate archr

#python HCA_ON/scripts/20_hotspot/5GO_bar_comb.py

R < HCA_ON/scripts/20_hotspot/5GO_bar_comb_update.R --no-save > HCA_ON/scripts/20_hotspot/5GO_bar_comb_update.out 2> HCA_ON/scripts/20_hotspot/5GO_bar_comb_update.err
