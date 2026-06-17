#!/bin/sh
source ~/.condainit
conda activate scarches
main="/dfs3b/ruic20_lab/junw42"
reference=${main}/HCA_ON/data/ref/scarches/retina
mkdir -p ${reference}
infile=${main}/HCA_ON/data/scvi/RNA/scvi/ONONH_hvg2k_epoch20_new_unknown_for_retina_raw.h5ad
outdir=${main}/HCA_ON/data/scvi/RNA/scvi/

# Download the reference model if not generated
#gdown -O “$reference/” https://zenodo.org/records/14014720/files/model.pt

# Download query .h5ad file
#gdown https://zenodo.org/records/14014720/files/query.h5ad

scarches2queryannobyscanvi -d "$outdir" -b query -e scarches -r "$reference" -k sampleid  -- "$infile"
