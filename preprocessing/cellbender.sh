#!/bin/bash
# https://cellbender.readthedocs.io/en/latest/reference/index.html


## ADD THIS PARAMETER WHEN ANALYSING ATAC DATA
##--projected-ambient-count-threshold 2 \ 
## FOR THE NEXT MULTIOME ANALYSIS FDR TO 0.05? maybe stick to 0.01

mamba activate cellbender

run_cellbender() {
    input_directory=$1
    for input_file in "$input_directory"/*/*/outs/raw_feature_bc_matrix.h5; do
        echo -e "Processing: $input_file"
        output_directory="${input_file%raw_feature_bc_matrix.h5}cellbender"
        mkdir -p "$output_directory"
        cd "$output_directory" || exit 1
        output_file="${input_file%raw_feature_bc_matrix.h5}cellbender/cellbender_feature_bc_matrix.h5"
        log_file="${input_file%raw_feature_bc_matrix.h5}cellbender/cellbender.log"
        echo -e "Output: $output_file\n\n"
        # running cellbender
        cellbender \
        remove-background \
            --cuda \
            --input "$input_file" \
            --output "$output_file" \
            --epochs 150 \
            --fpr 0.01 \
            > "$log_file" 2>&1

        exit_code=$?
        echo -e "\nexit code $exit_code" >> "$log_file"
        timestamp=$(date +"%D %T")
        echo "$timestamp" >> "$log_file"

    done
}

input_directory="/data/hadjantalab/lucas/sonja_project/preprocessing/mapping/scrnaseq925"
run_cellbender $input_directory