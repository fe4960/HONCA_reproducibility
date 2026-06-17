#!/bin/sh

#ls eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/*35UP.10DOWNMAGMA*ONONH_all_raw_normcount_only_5k_simple* | grep eye_QTL | sed -e "s/\://g"  > eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/MAGMA_list_major

ls eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/*35UP.10DOWNMAGMA*${1}* | grep eye_QTL | sed -e "s/\://g"  > eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/MAGMA_list_${1}


#ls eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/*35UP.10DOWNMAGMA*ONONH_all_raw_normcount_only_2k_celltype_rmUnk_simple* | grep eye_QTL | sed -e "s/\://g"  > eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/MAGMA_list_celltype


#ls eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/*35UP.10DOWNMAGMA*new5* | grep eye_QTL | sed -e "s/\://g"  > eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/MAGMA_list_Astrocyte
#ls eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/*35UP.10DOWNMAGMA_Oligodendrocyte* | grep eye_QTL | sed -e "s/\://g"  > eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/MAGMA_list_Oligodendrocyte
