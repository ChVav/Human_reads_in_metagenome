#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <DIR_IN> <DIR_OUT>"
    exit 1
fi

# Input and output directories from command line arguments
DIR_IN="$1"
DIR_OUT="$2"

# Create output directory if it doesn't exist
mkdir -p "$DIR_OUT"

# Exit code for tracking errors
exit_code=0

# Loop through all .fastq.gz files in the input directory
for file in "$DIR_IN"/*_1.fastq.gz; do

    # Extract the base name (removes directory path)
    basename=$(basename "$file")

    # Extract the sample ID (removes _1.fastq.gz)
    ID="${basename%_1.fastq.gz}"

    # Define the second fastq file based on the ID
    second_file="$DIR_IN/${ID}_2.fastq.gz"

    # Check if the second fastq file exists
    if [ ! -f "$second_file" ]; then
        echo "- Paired file for \"$second_file\" does not exist"
        exit_code=1
        break
    fi

    echo "- Found files: $file and $second_file"

    # Define output files for fastp
    first_file_out="$DIR_OUT/${ID}_fastp_1.fastq.gz"
    second_file_out="$DIR_OUT/${ID}_fastp_2.fastq.gz"
    failed_file="$DIR_OUT/${ID}_failed_fastp"

    echo "- Preprocessing and quality checking with fastp"

    # Run fastp
    singularity exec ./cleanreads.sif fastp \
    -i "$file" \
    -I "$second_file" \
    -o "$first_file_out" \
    -O "$second_file_out" \
    --failed_out "$failed_file" \
    --trim_poly_g \
    --detect_adapter_for_pe \
    --adapter_fasta ./adapters.fa \
    --cut_front 1 \
    --cut_tail 1 \
    --length_required 60 \
    -w 16 \
    -R "$ID"

    # Define output files for BBMap
    first_final_out="$DIR_OUT/${ID}_clean1.fastq.gz"
    second_final_out="$DIR_OUT/${ID}_clean2.fastq.gz"
    first_human_out="$DIR_OUT/${ID}_human1.fastq.gz"
    second_human_out="$DIR_OUT/${ID}_human2.fastq.gz"

    echo "- Removing human reads with BBMap"

    # Run BBMap
    singularity exec ./cleanreads.sif bbmap.sh t=16 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim in1="$first_file_out" in2="$second_file_out" outu1="$first_final_out" outu2="$second_final_out" outm1="$first_human_out" outm2="$second_human_out" 

done

# Exit with the final exit code
exit $exit_code
