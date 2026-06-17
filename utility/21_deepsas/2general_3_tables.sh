#!/bin/bash

dir=/dfs3b/ruic20_lab/junw42/software/deepsas/

cd $dir
source .venv/bin/activate

ct=$1
end=$2
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas/"
# Define input and output directories
SUBSAMPLE_DIR="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas/subsamples/" #"/bmbl_data/ahmed/eye_atlas/subsamples"  # Update this to your subsample directory
OUTPUT_DIR="${dir1}/outputs/" #"outputs"

#OUTPUT_DIR="/dfs3b/ruic20_lab/junw42/software/deepsas/outputs/" #"outputs"
DEVICE_INDEX=0  # Set your GPU device index
mkdir ${OUTPUT_DIR}
# Loop over all subsample files dynamically


for i in 0 #{0..3}
do
    # Define the file path dynamically
    FILE_PATH="${SUBSAMPLE_DIR}/subsample_${i}.h5ad"
#    FILE_PATH = "/dfs3b/ruic20_lab/junw42/software/deepsas/outputs/$"
    
    # Define the experiment name dynamically
    EXP_NAME="subsample_${i}"
    
    # Run your command (replace `your_command` with the actual command you want to execute)
    uv run python -u generate_3tables.py --input_data_count "$FILE_PATH" --output_dir $OUTPUT_DIR --exp_name $EXP_NAME --device_index 0 --retrain

    # Print (optional, for debugging)
    echo "Processing file: $FILE_PATH with experiment name: $EXP_NAME"
done





