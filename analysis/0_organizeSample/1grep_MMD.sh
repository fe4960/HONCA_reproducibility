#!/bin/sh
module load anaconda/2024.06

python HCA/ON_ONH/scripts/organize_sample/1grep_MMD.py

#for fl in ON_ONH_ATAC  ON_ONH
#do
#grep "MMD" HCA/ON_ONH/data/sample_list/${fl}.txt > HCA/ON_ONH/data/sample_list/ON_ONH_ATAC_MMD.txt
