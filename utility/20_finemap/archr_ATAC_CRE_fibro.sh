#!/bin/sh

source ~/.condainit

conda activate archr

R < HCA_ON/scripts/20_finemap/archr_ATAC_CRE_fibro.R --no-save > HCA_ON/scripts/20_finemap/archr_ATAC_CRE_fibro.out 2> HCA_ON/scripts/20_finemap/archr_ATAC_CRE_fibro.err
