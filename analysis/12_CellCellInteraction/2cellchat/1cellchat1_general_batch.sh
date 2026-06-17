#!/bin/sh

#######comd="HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general_cellphoneDB"


comd="HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general"
mkdir $comd
ct="celltype"
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/"

#name="Astro_oligo_opc_clean"
name="Astro_oligo_opc_clean_cellphoneDB"

wd="${main}/${name}"
mkdir $wd
cc=${wd}/celltype_list
#echo "Astrocyte" > $cc
#echo "Oligodendrocyte" >> $cc
#echo "Oligodendrocyte_precursor_cell" >> $cc
echo "Astro_ON" > $cc
echo "Astro_ONONH" >> $cc
echo "Astro_ONH" >> $cc
echo "Astro_retina" >> $cc
echo "OLIGO_SVEP1-hi" >> $cc
echo "OLIGO_LRRC7-hi" >> $cc
echo "OLIGO_RBFOX1-hi" >> $cc
echo "OLIGO_NAV2-hi" >> $cc
echo "OPC" >> $cc
#echo "OPC_RNF220+" >> $cc
#echo "OPC_GPR17+" >> $cc
#echo "OPC_LAMA2+" >> $cc

#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-3 -e $comd -j $name -- $wd $cc $ct

#name="Astro_endo_clean"
name="Astro_endo_clean_cellphoneDB"
#name="Astro_endo_mural"
wd="${main}/${name}"
mkdir $wd
cc=${wd}/celltype_list
#echo "Astrocyte" > $cc
#echo "Endothelial_cell" >> $cc
#echo "Mural_cell" >> $cc
echo "Astro_ON" > $cc
echo "Astro_ONONH" >> $cc
echo "Astro_ONH" >> $cc
echo "Astro_retina" >> $cc
echo "Capillary" >> $cc
echo "Vein/Venule" >> $cc
echo "Artery/Arteriole" >> $cc
echo "Venule_postcapillary" >> $cc
echo "Capillary_fenestrated" >> $cc
echo "Venule_POSTN+" >> $cc

#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-3 -e $comd -j $name -- $wd $cc $ct

#name="Astro_fibro_clean"

name="Astro_fibro_clean_cellphoneDB"

wd="${main}/${name}"
mkdir $wd
cc=${wd}/celltype_list
#echo "Astrocyte" > $cc
#echo "Fibroblast" >> $cc

echo "Astro_ON" > $cc
echo "Astro_ONONH" >> $cc
echo "Astro_ONH" >> $cc
echo "Astro_retina" >> $cc
echo "Fibro_arachnoid" >> $cc
echo "Fibro_pia" >> $cc
echo "Fibro_dura" >> $cc
echo "Fibro_RPEchoroid" >> $cc
echo "Fibro_lamina_cribrosa" >> $cc
echo "Fibro_sclera" >> $cc
echo "Fibro_perivascular" >> $cc

#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-3 -e $comd -j $name -- $wd $cc $ct



#sbatch -p free-gpu --mem=50GB --time=0-3 --error=${wd}/log.err --output=${wd}/log.out HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.sh$wd $cc $ct
#ct="subclass"
#name="ONH_age_12_cellphoneDB"

ct="subclass1"
#name="ONH_age_12_cellphoneDB_subclass1"
#name="ONH_age_12_subclass1"
#name="ONH_age_12_subclass1_fibro"

#name="ONH_age_12_subclass1_fibro_RPE"

name="ONH_age_12"
wd="${main}/${name}"
mkdir $wd
cc=${wd}/celltype_list

echo "Astro_ONHON" > $cc 

#echo "Astro_ONHON1_SLC4A11+MARCH1+" > $cc     
#echo "OLIGO1_SVEP1+"   >>  $cc
#echo "Fibroblast"      >> $cc
echo "Astro_ONH1_GABBR2+SV2B+" >> $cc
echo "Astro_ONH2_PAX5+GABBR2+"  >>  $cc
#echo "Microglia" >>  $cc                         
echo "Endothelial_cell"       >>  $cc           
echo "Oligodendrocyte_precursor_cell" >>  $cc    
echo "OLIGO2_LRRC7+"       >>  $cc              
echo "Mural_cell"          >>  $cc
#echo "OLIGO1"              >>  $cc               
#echo "Astro_ONHON2_CST3+APOE+"   >>  $cc         
#echo "Astro_retina1_PAX5+ME1+"   >>  $cc         
#echo "OLIGO1_RBFOX1+"           >>  $cc 

echo "Astro_retina"   >>  $cc  
echo "Astro_ON" >> $cc
echo "OLIGO1" >> $cc
echo "Fibro_sclera" >> $cc
echo "Fibro_pia" >> $cc
echo "Fibro_dura" >> $cc
echo "Fibro_RPEchoroid_SMOC2+" >> $cc
echo "Fibro_RPEchoroid_BMP5+" >> $cc
echo "Fibro_perivascular" >> $cc
echo "RPE" >> $cc

main1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/"
rds=${main1}/ONONH_all_rawcount_subset_age0to12.rds
sh sb.sh -n hpc3-l18-01 -c ${comd}.sh -m 20 -p s -t 0-2 -e $comd -j $name -- $wd $cc $ct $rds


#name="ONH_age_20to40_cellphoneDB"
name="ONH_age_20to40"

#name="ONH_age_20to40"
wd="${main}/${name}"
mkdir $wd
cc=${wd}/celltype_list

echo "Astro_ONHON1_SLC4A11+MARCH1+" > $cc     
echo "OLIGO1_SVEP1+"   >>  $cc
echo "Fibroblast"      >> $cc
echo "Astro_ONH1_GABBR2+SV2B+" >> $cc
echo "Astro_ONH2_PAX5+GABBR2+"  >>  $cc
echo "Microglia" >>  $cc                         
echo "Endothelial_cell"       >>  $cc           
echo "Oligodendrocyte_precursor_cell" >>  $cc    
echo "OLIGO2_LRRC7+"       >>  $cc              
echo "Mural_cell"          >>  $cc               
echo "OLIGO1"              >>  $cc               
echo "Astro_ONHON2_CST3+APOE+"   >>  $cc         
echo "Astro_retina1_PAX5+ME1+"   >>  $cc         
echo "OLIGO1_RBFOX1+"           >>  $cc 

main1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/"
rds=${main1}/ONONH_all_rawcount_subset_age20to40.rds
sh sb.sh -n hpc3-l18-00 -c ${comd}.sh -m 20 -p s -t 0-2 -e $comd -j $name -- $wd $cc $ct $rds

#name="ONH_age_70_cellphoneDB"

name="ONH_age_70"
wd="${main}/${name}"
mkdir $wd
cc=${wd}/celltype_list

echo "Astro_ONHON1_SLC4A11+MARCH1+" > $cc     
echo "OLIGO1_SVEP1+"   >>  $cc
echo "Fibroblast"      >> $cc
echo "Astro_ONH1_GABBR2+SV2B+" >> $cc
echo "Astro_ONH2_PAX5+GABBR2+"  >>  $cc
echo "Microglia" >>  $cc                         
echo "Endothelial_cell"       >>  $cc           
echo "Oligodendrocyte_precursor_cell" >>  $cc    
echo "OLIGO2_LRRC7+"       >>  $cc              
echo "Mural_cell"          >>  $cc               
echo "OLIGO1"              >>  $cc               
echo "Astro_ONHON2_CST3+APOE+"   >>  $cc         
echo "Astro_retina1_PAX5+ME1+"   >>  $cc         
echo "OLIGO1_RBFOX1+"           >>  $cc 

main1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/"
rds=${main1}/ONONH_all_rawcount_subset_age70.rds
sh sb.sh -n hpc3-l18-00 -c ${comd}.sh -m 40 -p s -t 0-2 -e $comd -j $name -- $wd $cc $ct $rds
