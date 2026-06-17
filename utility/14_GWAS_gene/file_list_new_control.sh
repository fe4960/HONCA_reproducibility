ls  HCA_OCS/data/2_GWAS/LDSC_hg38/ldsc_result/ | grep cell_type | awk -F "." '{print $1}'  > HCA_ON/scripts/14_GWAS_gene//file_list_new_control
