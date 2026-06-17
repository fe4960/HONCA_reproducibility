#!/bin/sh
module load R/4.4.2 
R < HCA_ON/scripts/7_DEG/11eQTL_RNASEH2B_RGC.R --no-save > HCA_ON/scripts/7_DEG/11eQTL_RNASEH2B_RGC.out 2> HCA_ON/scripts/7_DEG/11eQTL_RNASEH2B_RGC.err
