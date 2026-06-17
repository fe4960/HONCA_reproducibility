#!/bin/sh
comd="HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat3_general"
mkdir $comd
ct="celltype"
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11cellchat/"

name="Astro_oligo_opc_clean"
wd="${main}/${name}"
mkdir $wd
lr=${wd}/lr_list

#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-3 -e $comd -j $name -- $wd $cc $ct

name="Astro_endo_clean"

#name="Astro_endo_mural"
wd="${main}/${name}"
mkdir $wd
cc=${wd}/pathway

sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-3 -e $comd -j $name -- $wd $cc

name="Astro_fibro_clean"
wd="${main}/${name}"
mkdir $wd
cc=${wd}/pathway

#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-3 -e $comd -j $name -- $wd $cc $ct



#sbatch -p free-gpu --mem=50GB --time=0-3 --error=${wd}/log.err --output=${wd}/log.out HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.sh$wd $cc $ct

