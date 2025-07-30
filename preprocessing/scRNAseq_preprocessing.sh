## QC FASTQ -------------------------------------------------------
input_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/raw/scrnaseq925"
output_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/fastqc/scrnaseq925"

for sample_path in "$input_dir"/*; do
    [ -d "$sample_path" ] || continue  # Skip non-directories
    sample=$(basename "$sample_path")
    mkdir -p "$output_dir/$sample"
    find "$sample_path" -name "*.fastq.gz" | \
        parallel --jobs 4 fastqc --outdir "$output_dir/$sample" {}
done

## QC MULTIQC all files
for sample in $(ls $output_dir)
do
    multiqc --outdir $output_dir/$sample/multiqc_all_reads $output_dir/$sample/*
done


## QC MULTIQC biological reads only
for sample in $(ls $output_dir)
do
    multiqc --outdir $output_dir/$sample/multiqc_bio_reads $output_dir/$sample/*_R2_*
done


### MAPPING CELLRANGER PARALLEL (WITH libraries.csv) ---------------------------------

# Function to run cellranger count for one sample
process_sample() {
    sample=$1
    libraries_dir=$2
    output_dir=$3
    transcriptome=$4
    num_cores=$5
    localmem=$6

    library="$libraries_dir/$sample.csv"
    sample_output_path="$output_dir/$sample"
    log_file="$sample_output_path/$sample.log"

    mkdir -p "$sample_output_path"
    cd "$sample_output_path" || exit 1

    /home/schroel1/cellranger-8.0.1/cellranger count \
        --id="$sample" \
        --nosecondary \
        --create-bam=true \
        --transcriptome="$transcriptome" \
        --libraries="$library" \
        --localcores="$num_cores" \
        --localmem="$localmem" \
        > "$log_file" 2>&1

    exit_code=$?
    echo -e "\nexit code $exit_code" >> "$log_file"
    timestamp=$(date +"%D %T")
    echo "$timestamp" >> "$log_file"
    if [ $exit_code -ne 0 ]; then
        echo "Error: Cellranger failed for sample $sample" >> "$log_file"
    fi
}

export -f process_sample

# Define paths
libraries_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/libraries/scrnaseq925"
output_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/mapping/scrnaseq925"
transcriptome="/admin/opt/common/CentOS_7/cellranger/cellranger-7.1.0/refdata-cellranger/refdata-gex-mm10-2020-A"

# Compute resources
num_cores=$(( $(nproc) / 3 ))   
localmem=$(( 300 / 3 ))         

# Run in parallel over all .csv files in the libraries_dir
find "$libraries_dir" -name "*.csv" -exec basename {} .csv \; | \
parallel --jobs 3 process_sample {} "$libraries_dir" "$output_dir" "$transcriptome" "$num_cores" "$localmem"



