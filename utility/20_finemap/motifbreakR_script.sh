#!/bin/sh

module load R/4.4.2
R < HCA_ON/scripts/20_finemap/motifbreakR_script.R --no-save > HCA_ON/scripts/20_finemap/motifbreakR_script.out 2> HCA_ON/scripts/20_finemap/motifbreakR_script.err 

