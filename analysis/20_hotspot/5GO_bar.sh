#!/bin/sh

source ~/.condainit

conda activate qtl


#python HCA_ON/scripts/20_hotspot/5GO_bar.py

R < HCA_ON/scripts/20_hotspot/5GO_bar.R --no-save > HCA_ON/scripts/20_hotspot/5GO_bar.out 2> HCA_ON/scripts/20_hotspot/5GO_bar.err
