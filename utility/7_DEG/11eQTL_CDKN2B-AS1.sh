#!/bin/sh
module load R/4.4.2 
R < HCA_ON/scripts/7_DEG/11eQTL_CDKN2B-AS1.R --no-save > HCA_ON/scripts/7_DEG/11eQTL_CDKN2B-AS1.out 2> HCA_ON/scripts/7_DEG/11eQTL_CDKN2B-AS1.err
