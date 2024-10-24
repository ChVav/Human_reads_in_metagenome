#!/bin/bash
# With human reads already extracted from metagenomes, can I still recover XY chromosome info?

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

for fastq1 in "$DIR_IN"/*_human1.fastq.gz; do

    # Extract the base name (removes directory path)
    basename=$(basename "$fastq1")

    # Extract the sample ID (removes _human1.fastq.gz)
    ID="${basename%_human1.fastq.gz}"

    # Define the second fastq file based on the ID
    fastq2="$DIR_IN/${ID}_human2.fastq.gz"

    # Check if the second fastq file exists
    if [ ! -f "$fastq2" ]; then
        echo "Paired file for $fastq1 not found!"
        continue
    fi

    # Define the output files and paths
    output_sam="$DIR_OUT/${ID}.sam"
    output_bam="$DIR_OUT/${ID}.bam"
    output_bam_index="$DIR_OUT/${ID}.bam.bai"
    output_bam2="$DIR_OUT/${ID}_rg.bam"
    output_bam2_index="$DIR_OUT/${ID}_rg.bai"
    output_bam3="$DIR_OUT/${ID}_rg_nodup.bam"
    metrics="$DIR_OUT/${ID}_rg_nodup.metrics.txt"
    out_fastq1="$DIR_OUT/rm_dup_${ID}_mapped_R1.fastq"
    out_fastq2="$DIR_OUT/rm_dup_${ID}_mapped_R2.fastq"
    xynonpar_summary="$DIR_OUT/bedcov_${ID}_nonPAR_XY.txt"

    # Run bowtie2 with the extracted ID
    singularity exec bowtie2_samtools.sif bowtie2 --threads 16 -R 3 --no-discordant \
        -x human_g1k_v37 \
        -1 "$fastq1" \
        -2 "$fastq2" \
        -S "$output_sam"

    echo "Bowtie2 alignment completed for $ID, output saved to $output_sam"

    # Extract only properly mapped reads to bam and index
    singularity exec bowtie2_samtools.sif samtools view -Sb -@ 2 "$output_sam" | samtools sort -@ 2 -m 20G | \
    samtools view -u -f 3 > "$output_bam"
    rm "$output_sam"
    singularity exec bowtie2_samtools.sif samtools index "$output_bam"

    # Add read groups
    singularity exec ./picard.sif java -Xmx16g -jar /usr/local/bin/picard/picard.jar \
    AddOrReplaceReadGroups I="$output_bam" \
    O="$output_bam2" \
    RGID=FLOWCELLID \
    RGLB=${ID}_test \
    RGPU=DUMMY \
    RGPL=illumina \
    RGSM=${ID} \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

    # Remove duplicates
    singularity exec ./picard.sif java -Xmx16g -jar /usr/local/bin/picard/picard.jar \
    MarkDuplicates I="$output_bam2" \
    O="$output_bam3" \
    M="$metrics" \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT

    # Index the final BAM file
    singularity exec bowtie2_samtools.sif samtools index "$output_bam3"

    # Extract the non-duplicate reads
    singularity exec ./picard.sif java -Xmx16g -jar /usr/local/bin/picard/picard.jar \
    SamToFastq I="$output_bam3" \
    F="$out_fastq1" F2="$out_fastq2"

    # Clean up
    gzip -f "$out_fastq1"
    gzip -f "$out_fastq2"
    rm "$output_bam" "$output_bam_index" "$output_bam2" "$output_bam2_index"

    # Check XY non-par regions
    singularity exec bowtie2_samtools.sif samtools view -u \
    "$output_bam3" | \
    singularity exec bedtools.sif bedtools coverage -a nonpar.bed \
    -b stdin -counts -F 0.1 | \
    awk 'BEGIN{OFS="\t"}{print $0"\t"$4/($3-$2)"\t"$4}' | \
    sed '1i CHR\tSTART\tEND\tBP\tDepth\tCount' > "$xynonpar_summary"

done
