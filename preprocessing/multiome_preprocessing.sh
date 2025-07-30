
## some renaming stuff 
prefix="Fauci2_0040_" 
for file in *; do mv "$file" "${prefix}${file}"
done

postfix="Bono_0004"
for file in *; do
    if [ -f "$file" ]; then
        base="${file%%.*}"
        ext="${file#"$base"}"
        mv "$file" "${base}_${postfix}${ext}"
    fi
done

postfix="_Fauci2_0040.fastq.gz"
prefix="Fauci2_0040_"
for file in *$postfix; do
    if [ -f "$file" ]; then
        base="${file%$postfix}"  # Remove the _Bono_0007.fastq.gz part
       mv "$file" "${prefix}${base}.fastq.gz"
    fi
done

for file in *; do
    if [[ -f "$file" && "$file" == *-* ]]; then
        newname="${file//-/_}"
        echo "$file" "$newname"
    fi
done


pattern="_Bono_0004"
for file in *; do
    if [[ -f "$file" && "$file" ==  *$pattern* ]]; then
        newname="${file//$pattern/}"
        mv "$file" "$newname"
    fi
done

## the fucking naming of these files is giving me a headache because they are run on two sequencers,
# they are run on either bono_0007 or fauci2_0040 
## I can still dsitinguish them via rest of the nomenclature tho, so I will remove the info of the machine.
## in the raw data the structure of the raw files is still preserved so It is obvious which reads come from which sequencer


## QC FASTQ -------------------------------------------------------
input_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/raw/multiome95/"
output_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/fastqc/multiome95"
for sample in $(ls $input_dir)
do  
    mkdir --parents $output_dir/$sample
    ls $input_dir/$sample/*.fastq.gz | parallel -I % --jobs $(ls $input_dir/$sample/ | wc -l) fastqc --outdir $output_dir/$sample %
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


## MAPPING CELLRANGER-ARC SEQUENTIALLY

map_multiome_samples() {
    libraries_dir=$1
    for file in "$libraries_dir"/*.csv; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            sample="${filename%.csv}"           
            output="$2/$sample"
            library="$libraries_dir/$sample.csv"
            reference="$3"
            num_cores=$(( 20 ))
            localmem=$(( 300 ))
            log_file="$output/$sample.log"
            module load cellrangerARC/2.0.2
            mkdir -p "$output"
            cd "$output" || exit 1
            cellranger-arc count \
            --id="$sample" \
            --reference="$reference" \
            --libraries="$library" \
            --localcores="$num_cores" \
            --localmem="$localmem" \
            --uiport=3600 \
            > $log_file 2>&1

            exit_code=$?
            echo -e "\nexit code $exit_code" >> "$log_file"
            timestamp=$(date +"%D %T")
            echo "$timestamp" >> "$log_file"
            # Optional error message in log if the process fails
            if [ $exit_code -ne 0 ]; then
                echo "Error: Cellranger failed for sample $sample" >> "$log_file"
            fi
        fi
    done
}


libraries_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/libraries"
output_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/mapping"
reference="/data/hadjantalab/cellranger_refData/refdata-cellranger-arc-mm10-2020-A-2.0.0"

map_multiome_samples "$libraries_dir" "$output_dir" "$reference"

## MAPPING CELLRANGER PARALLEL

# Function to run cellranger-arc for one sample
process_sample() {
    sample=$1
    libraries_dir=$2
    output_dir=$3
    reference=$4

    library="$libraries_dir/$sample.csv"
    output="$output_dir/$sample"
    log_file="$output/$sample.log"
    num_cores=$5
    localmem=$6

    mkdir -p "$output"
    cd "$output" || exit 1

    cellranger-arc count \
        --id="$sample" \
        --reference="$reference" \
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
libraries_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/libraries/multiome95"
output_dir="/data/hadjantalab/lucas/sonja_project/preprocessing/mapping/multiome95"
reference="/data/hadjantalab/cellranger_refData/refdata-cellranger-arc-mm10-2020-A-2.0.0"
num_cores=$(( (60) / 3 ))
localmem=$(( 300 / 3 ))

# Run in parallel over all CSV filenames
find "$libraries_dir" -name "*.csv" -exec basename {} .csv \; | \
parallel --jobs 3 process_sample {} "$libraries_dir" "$output_dir" "$reference" "$num_cores" "$localmem"



