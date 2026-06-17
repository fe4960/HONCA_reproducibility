#!/bin/bash

PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:128

dir=/dfs3b/ruic20_lab/junw42/software/deepsas/

cd $dir
source .venv/bin/activate

ct=$1
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas/"
# Define input and output directories
SUBSAMPLE_DIR="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas/subsamples/" #"/bmbl_data/ahmed/eye_atlas/subsamples"  # Update this to your subsample directory
OUTPUT_DIR="${dir1}/outputs"
#cd $dir1
DEVICE_INDEX=0  # Set your GPU device index
#mkdir ${OUTPUT_DIR}
# Loop over all subsample files dynamically
for FILE_PATH in "$SUBSAMPLE_DIR"/subsample_*.h5ad; do
    # Extract the iteration number from the filename (e.g., "subsample_1.h5ad" -> "1")
    FILE_NAME=$(basename "$FILE_PATH")  # Extract filename
#    ITER_NUM=$(echo "$FILE_NAME" | grep -oP '\d+')  # Extract numeric part
    EXP_NAME=$(echo "$FILE_NAME" | cut -d "." -f 1)	
    # Define the experiment name dynamically
#    EXP_NAME="Data_${ITER_NUM}"


    # Run the DeepSAS command for each subsample
    uv run python -u deepsas_v1.py \
        --input_data_count "$FILE_PATH" \
        --output_dir "$OUTPUT_DIR" \
        --exp_name "$EXP_NAME" \
        --device_index "$DEVICE_INDEX" \
        --retrain \
	--n_genes 5000

#	--batch_id "sampleid" \

    # Print for debugging
    echo "Processing file: $FILE_PATH with experiment name: $EXP_NAME"
done




