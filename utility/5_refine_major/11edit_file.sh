#!/bin/sh

for i in $( ls HCA_ON/data/5_refine_major/scvi/major/figures/*hvg10000_nonsb_epochnone_seurat_v3_rs_1_clean_ON_ONH_new*png )
do
fn=$(echo $i | sed "s/nonsb/sb/g" )
echo $fn
mv $i $fn
done


for i in $( ls HCA_ON/data/5_refine_major/scvi/major/*hvg10000_nonsb_epochnone_seurat_v3_rs_1_clean_ON_ONH_new* )
do
fn=$(echo $i | sed "s/nonsb/sb/g" )
echo $fn
mv $i $fn
done
