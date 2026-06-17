#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

#err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution"
err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution_sb"

mkdir -p $err

#for cell in AC Astrocyte BC Fibroblast MG Macrophage Melanocyte Microglia Oligodendrocyte RGC RPE  
#for cell in Astrocyte
#for cell in Fibroblast
#for cell in  Astrocyte Cone HC MG Macrophage Mural_cell Schwann_cell  NK_T_cell  Oligodendrocyte_precursor_cell Microglia Endothelial_cell AC BC RGC RPE Fibroblast
for cell in major
#for cell in Oligodendrocyte
#for cell in Schwann_cell
#for cell in Melanocyte
#for cell in Macrophage
#for cell in Endothelial_cell 
#for cell in NK_T_cell
#for cell in Oligodendrocyte_precursor_cell
#for cell in RPE
#for cell in MG
#for cell in Mural_cell
do	
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean
indir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/

#for res in 0.1 0.2 0.4 #0.3 1  #0.4 0.5 0.6
#for res in 0.9 1  #0.3 0.4 0.5 0.6 
#for res in 0.7 0.8 0.9
#for res in  0.8 0.9
#for res in 0.2 0.3 0.5 0.6
for res in 0.4 0.5 0.6 0.7 0.8
do

#h5ad=${outdir}/${cell}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_rm_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm_rm_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg20000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_16_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg2000_epochnone_seurat_rs_0.4_clean_rm4_15_sb_scvi_trg.h5ad


h5ad=${outdir}/${cell}_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_scvi_trg.h5ad


#h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_v3_rs_1_span_1_subclass_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg2000_epochnone_seurat_v3_rs_0.4_clean_rm_rm_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg2000_epochnone_seurat_v3_rs_0.4_clean_rm_scvi_trg.h5ad

#bname=${cell}_res_${res}

#bname=${cell}_res_${res}_clean
#bname=${cell}_res_${res}_clean_rm
#bname=${cell}_res_${res}_clean_rm_rm
#bname=${cell}_res_${res}_clean_rm_rm_rm
#bname=${cell}_res_${res}_clean_rm4_15_16
#bname=${cell}_res_${res}_clean_sb
#bname=${cell}_res_${res}_clean_rm4_15
#bname=${cell}_res_${res}_clean_rm4_15_sb
bname=${cell}_hvg_2000_fl_${res}_res_${res}_ononh_retina_mg_sb
#outdir=$indir
rmlist="none" #"rm_cluster_rs" #"none" #rm_cluster_res_1
mk="none" #"sanes_mk"
save="n"
#save="t"
sbatch --mem=40GB -p free --account=ruic20_lab --time=0-2 --output=${err}/${cell}_${res}.out --error=${err}/${cell}_${res}.err  HCA_ON/scripts/5_refine_major/5refine_resolution.sh $h5ad $res ${batch_key} $bname $outdir $rmlist $mk $save $indir
done
done
