#!/bin/sh

module load bedtools2/2.30.0

#peak_file="/storage/chenlab/Users/junwang/sc_human_retina/data/snATAC_seq_peak/narrowPeaks_macs3.combined.bed.mainChr_rmBlacklist_rmY_flt_2fpkm"
#peak_file="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/peakset_clusters_major_full_peak_anno_hg38_hg19"

#sed -e "s/\:/\t/g" $peak_file | sed -e "s/\-/\t/g" |  sed -e "s/chr//g" > ${peak_file}_bed

#peak_file="/storage/chenlab/Users/junwang/human_ret_anc/data/proj1_clean_remap_HRCA_rmMulti_clean_final/peakset_clusters_major_full_peak_anno_hg38_hg19"

peak_file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/peakset_clusters_major_full_peak_anno_hg38_hg19"

sed -e "s/\:/\t/g" $peak_file | sed -e "s/\-/\t/g" |  sed -e "s/chr//g" | awk '{print $1"\t"$2-1"\t"$3}' > ${peak_file}_bed



#for i in clinical_c_Block_H30-H36 clinical_c_H33 clinical_c_H35 clinical_c_H36 selfReport
for i in  Meta_CA_age_sex_DA_caucasians_MAF0.01_20160807_newGC Meta_DA_age_sex_caucasians_20160806 Meta_IOP_age_sex_caucasians_MAF0.01_20160806 Meta_VCDR_age_sex_caucasians_withoutSouthampton_MAF0.01_20160806 UKBBc.IOP.txt
do 

#file=/storage/chen/home/jw29/sc_human_retina/data/GWAS/SummaryStat/geneAtlas/results/${i}/imputed.allWhites_all_chr_annotated_5e8_rmChr
file=reference/GWAS/finemap_sumstat/${i}_5e8
echo "##fileformat=VCFv4.2" > ${file}_format_atlas_new.vcf 
echo "#CHROM	POS	ID	REF	ALT	QUAL	INFO" >> ${file}_format_atlas_new.vcf 
#awk '{print $8"\t"$9"\t"$1"\t"$7"\t"$2"\t"$4"\t"$6}' ${file} >> ${file}_format_atlas.vcf 
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8}' ${file} >> ${file}_format_atlas_new.vcf #3 1 2 4 5 6 8 

bedtools intersect -wo -a ${file}_format_atlas_new.vcf -b ${peak_file}_bed  > ${file}_format_atlas_peak_ON

done


#file=sc_human_retina/data/GWAS/SummaryStat/FritscheLG2016_26691988/FritscheLG2016_26691988NG/Fritsche-26691988.txt_reform_5e8
for i in Khor-27064256_reform UKBBc.IOP.txt
do	
file=reference/GWAS/finemap_sumstat/${i}_5e8
echo "##fileformat=VCFv4.2" > ${file}_format_atlas_new.vcf
echo "#CHROM    POS     ID      REF     ALT     QUAL    INFO" >> ${file}_format_atlas_new.vcf
awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$6}' ${file} >>  ${file}_format_atlas_new.vcf #1 2 3 4 5 8 6

bedtools intersect -wo -a ${file}_format_atlas_new.vcf -b ${peak_file}_bed > ${file}_format_atlas_peak_ON
done



file=reference/GWAS/finemap_sumstat/GCST90011770_buildGRCh37_5e8

echo "##fileformat=VCFv4.2" > ${file}_format_atlas_new.vcf
echo "#CHROM    POS     ID      REF     ALT     QUAL    INFO" >> ${file}_format_atlas_new.vcf

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8}' ${file} >> ${file}_format_atlas_new.vcf #3 1 2 4 5 6 8 

bedtools intersect -wo -a ${file}_format_atlas_new.vcf -b ${peak_file}_bed > ${file}_format_atlas_peak_ON




