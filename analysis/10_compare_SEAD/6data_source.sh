#!/bin/sh
module load R/4.4.2

R < HCA_ON/scripts/10_compare_SEAD/6data_source.R --no-save > HCA_ON/scripts/10_compare_SEAD/6data_source.err 2> HCA_ON/scripts/10_compare_SEAD/6data_source.out 
