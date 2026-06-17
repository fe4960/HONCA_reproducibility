#!/bin/sh

source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/12ispatial_rest5.R --no-save > HCA_ON/scripts/17_xenium/12ispatial_rest5.out 2 > HCA_ON/scripts/17_xenium/12ispatial_rest5.err 
