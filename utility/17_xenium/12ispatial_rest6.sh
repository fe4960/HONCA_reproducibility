#!/bin/sh

source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/12ispatial_rest6.R --no-save > HCA_ON/scripts/17_xenium/12ispatial_rest6.out 2 > HCA_ON/scripts/17_xenium/12ispatial_rest6.err 
