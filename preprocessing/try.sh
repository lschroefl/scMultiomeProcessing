run_velocyto() {
    sample=$1
    samples_dir=$2
    reference=$3
    num_cores=$4
    memory_per_thread=$5

    output="$samples_dir/$sample/$sample/outs/velocyto"
    log_file="$output/$sample.log"

    bcfile="$samples_dir/$sample/$sample/outs/cellbender/cellbender_feature_bc_matrix_cell_barcodes.csv"
    bam_file="$samples_dir/$sample/$sample/outs/possorted_genome_bam.bam"

    echo "Using bcfile: $bcfile"
    echo "Using bam_file: $bam_file"

    if [[ ! -f "$bcfile" ]]; then
        echo "Error: Barcode file does not exist: $bcfile"
        return 1
    fi

    mkdir -p "$output"
    cd "$output" || exit 1

    velocyto run -b "$bcfile" \
        --outputfolder "$output" \
        --samtools-threads "$num_cores" \
        --samtools-memory "$memory_per_thread" \
        --verbose \
        "$bam_file" \
        "$reference" > "$log_file" 2>&1

    exit_code=$?
    echo -e "\nexit code $exit_code" >> "$log_file"
    timestamp=$(date +"%D %T")
    echo "$timestamp" >> "$log_file"
    if [ $exit_code -ne 0 ]; then
        echo "Error: Failed for sample $sample" >> "$log_file"
    fi
}


export -f run_velocyto

# Define paths
libraries_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/libraries/scrnaseq925"
samples_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/mapping/scrnaseq925"
reference="/data/hadjantalab/cellranger_refData/mm10-2020-ref/genes/genes.gtf"
num_cores=$(( (60) / 3 ))
memory_per_thread=$(( (300 / 3) / num_cores ))

# Run in parallel over all CSV filenames
find "$libraries_dir" -name "*.csv" -maxdepth 1 -exec basename {} .csv \; | \
parallel --jobs 3 run_velocyto {} "$samples_dir" "$reference" "$num_cores" "$memory_per_thread"
