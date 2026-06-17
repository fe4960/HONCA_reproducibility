#!/bin/sh
main=/dfs3b/ruic20_lab/junw42
err=${main}/HCA_ON/scripts/14_GWAS_gene/2run_MAGMA
mkdir $err

#s=ONONH_snRNA_subclass_5k
#s=RPEchoroid_scRNA_majorclass_5k

#s=ONONH_snRNA_majorclass_5k
#s=retina_snRNA_majorclass_5k
#s=TMCB_snRNA_subtype_5k
#l=majorclass
#l=subclass
#l=subclass
#for s in RPEchoroid_scRNA_majorclass_5k RPEchoroid_snRNA_majorclass_5k
#for s in retina_snRNA_majorclass_5k TMCB_snRNA_subtype_5k
#for s in TMCB_snRNA_subtype_5k
#do
c=Astrocyte
s=Astrocyte_subclass_new4_clean_simple
#sbatch -p free-gpu --account=ruic20_lab --time=0-1 --mem=20GB --error=${err}/${s}.err --output=${err}/${s}.out ${err}.sh $s $l $c
l=subtype
#c=Fibro_CB
#s=Fibro_CB
c=Fibro_CB
s=cb_fibro_subtype_02_rawcount_anno_clean_simple
sbatch -p standard --account=ruic20_lab --time=0-2 --mem=20GB --error=${err}/${s}.err --output=${err}/${s}.out ${err}.sh $s $l $c

c=Fibro_TM
s=Fibro_TM
########sbatch -p free-gpu --account=ruic20_lab --time=0-1 --mem=20GB --error=${err}/${s}.err --output=${err}/${s}.out ${err}.sh $s $l $c


#done
#s=ONONH_snRNA_majorclass_5k
#l=majorclass

#sbatch -p standard --account=ruic20_lab --time=0-6 --mem=20GB --error=${err}/${s}.err --output=${err}/${s}.out ${main}/eye_QTL/scripts/2_GWAS/MAGMA/2run_MAGMA.sh $s $l
