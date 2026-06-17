#!/bin/sh

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/"
#export PERL5LIB=/storage/chen/home/jw29/software/perl/lib/perl5/Math-CDF-0.1/lib/perl5/x86_64-linux/:$PERL5LIB

err=${dir}/9anno_pip_perl
mkdir $err

main="reference/GWAS/finemap_sumstat/"

file="/dfs3b/ruic20_lab/junw42/${main}/GCST90011770_buildGRCh37_5e8"
n="Gharahkhani_33627673NC"
sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/${n}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno $n





file="/dfs3b/ruic20_lab/junw42/${main}/Khor-27064256_reform_5e8"
n="Khor2016_27064256PCAG"
#sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/POAG.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno Khor-27064256

sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/${n}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno $n

file="/dfs3b/ruic20_lab/junw42/${main}/Meta_CA_age_sex_DA_caucasians_MAF0.01_20160807_newGC_5e8"
n="Springelkamp2017_28073927CA"
#sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/POAG.out --error=${err}/POAG.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno Meta_CA
sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/${n}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno $n

file="/dfs3b/ruic20_lab/junw42/${main}/Meta_DA_age_sex_caucasians_20160806_5e8"
n="Springelkamp2017_28073927DA"
#sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/POAG.out --error=${err}/POAG.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno Meta_DA
sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/${n}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno $n

file="/dfs3b/ruic20_lab/junw42/${main}/Meta_IOP_age_sex_caucasians_MAF0.01_5e8"
n="Springelkamp2017_28073927IOP"
#########sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/POAG.out --error=${err}/POAG.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno Meta_IOP
sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/${n}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno $n

file="/dfs3b/ruic20_lab/junw42/${main}/Meta_VCDR_age_sex_caucasians_withoutSouthampton_MAF0.01_20160806_5e8"
n="Springelkamp2017_28073927VCDR"
##########sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/POAG.out --error=${err}/POAG.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno Meta_VCDR
sbatch --account=ruic20_lab   -p standard --mem=10GB --output=${err}/${n}.out --error=${err}/${n}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_finemap/9anno_pip_perl.sh  ${file}_format_atlas_peak_ON  ${file}_anno $n


