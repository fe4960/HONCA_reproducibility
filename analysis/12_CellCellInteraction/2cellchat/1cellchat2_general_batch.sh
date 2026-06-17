#!/bin/sh
comd="HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat2_general" #_rest"
mkdir $comd
ct="celltype"
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/"

#name="Astro_oligo_opc"
name="Astro_oligo_opc_clean"

wd="${main}/${name}"
mkdir $wd

#########sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-4 -e $comd -j $name -- $wd 

#name="Astro_endo_mural"
#name="Astro_endo_clean"
#name="ONH_age_12"
#name="ONH_age_12_cellphoneDB"
#name="ONH_age_12_cellphoneDB_subclass1"
#name="ONH_age_12_subclass1"
#name="ONH_age_12_subclass1_fibro"
#name="ONH_age_12_subclass1_fibro_RPE"

name="ONH_age_12"
wd="${main}/${name}"
mkdir $wd

sh sb.sh -n hpc3-l18-01 -c ${comd}.sh -m 50 -p s -t 0-3 -e $comd -j $name -- $wd

#name="Astro_fibro"
#name="Astro_fibro_clean"
name="ONH_age_20to40"
#name="ONH_age_20to40_cellphoneDB"

wd="${main}/${name}"
mkdir $wd

sh sb.sh -n hpc3-l18-00 -c ${comd}.sh -m 50 -p s -t 0-3 -e $comd -j $name -- $wd

#name="ONH_age_70_cellphoneDB"

name="ONH_age_70"
wd="${main}/${name}"
mkdir $wd

sh sb.sh -n hpc3-l18-00 -c ${comd}.sh -m 50 -p s -t 0-3 -e $comd -j $name -- $wd



#sbatch -p free-gpu --mem=50GB --time=0-3 --error=${wd}/log.err --output=${wd}/log.out HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.sh$wd $cc $ct

