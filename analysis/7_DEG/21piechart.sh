#!/bin/sh
source ~/.condainit
conda activate qtl

R < HCA_ON/scripts/7_DEG/21piechart.R --no-save > HCA_ON/scripts/7_DEG/21piechart.out 2> HCA_ON/scripts/7_DEG/21piechart.err
