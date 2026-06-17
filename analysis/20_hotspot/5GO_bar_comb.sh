#!/bin/sh

source ~/.condainit

conda activate qtl


#python HCA_ON/scripts/20_hotspot/5GO_bar_comb.py

R < HCA_ON/scripts/20_hotspot/5GO_bar_comb1.R --no-save > HCA_ON/scripts/20_hotspot/5GO_bar_comb1.out 2> HCA_ON/scripts/20_hotspot/5GO_bar_comb1.err

#R < HCA_ON/scripts/20_hotspot/5GO_bar_comb.R --no-save > HCA_ON/scripts/20_hotspot/5GO_bar_comb.out 2> HCA_ON/scripts/20_hotspot/5GO_bar_comb.err
