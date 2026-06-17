#!/bin/sh

module load R/4.4.2
R < eye_QTL/scripts/2_GWAS/MAGMA/5dotplot_clean.R --no-save > eye_QTL/scripts/2_GWAS/MAGMA/5dotplot_clean.out 2> eye_QTL/scripts/2_GWAS/MAGMA/5dotplot_clean.err
