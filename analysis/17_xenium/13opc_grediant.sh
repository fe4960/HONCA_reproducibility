#!/bin/sh
module load R/4.4.2 

R < HCA_ON/scripts/17_xenium/13opc_grediant.R --no-save > HCA_ON/scripts/17_xenium/13opc_grediant.out 2> HCA_ON/scripts/17_xenium/13opc_grediant.err
