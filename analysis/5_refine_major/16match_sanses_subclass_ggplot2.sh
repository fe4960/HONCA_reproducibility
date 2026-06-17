#!/bin/sh

module load R/4.4.2

R --vanilla --args  $1 $2 < HCA_ON/scripts/5_refine_major/16match_sanses_subclass_ggplot2.R  
