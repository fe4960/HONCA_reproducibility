#!/bin/sh

source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/12ispatial.R --no-save > HCA_ON/scripts/17_xenium/12ispatial.out 2 > HCA_ON/scripts/17_xenium/12ispatial.err 
