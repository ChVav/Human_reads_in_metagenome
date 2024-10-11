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

    # Define the output SAM file path
    output_sam="$DIR_OUT/${ID}.sam"

    # Run bowtie2 with the extracted ID
    singularity exec bowtie2_samtools.sif bowtie2 -R 3 --no-discordant \
        -x human_g1k_v37.fasta.gz \
        -1 "$fastq1" \
        -2 "$fastq2" \
        -S "$output_sam"

    echo "Bowtie2 alignment completed for $ID, output saved to $output_sam"

done

    # Extract only properly mapped reads
    #samtools view -Sb -@ 2 ${ID}.sam | samtools sort -@ 2 -m 20G | \
    #samtools view -u -f 3 > ${ID}_mapped_pre.bam
    #rm ${ID}.sam
    #samtools index ${ID}_mapped_pre.bam

    # Add read groups
    #java -Xmx16g -jar /path/to/picard/picard.jar \
    #AddOrReplaceReadGroups I=${ID}_mapped_pre.bam \
    #O=${ID}_mapped.bam \
    #RGID=FLOWCELLID \
    #RGLB=${ID}_test \
    #RGPU=DUMMY \
    #RGPL=illumina \
    #RGSM=${ID} \
    #VALIDATION_STRINGENCY=LENIENT \
    #CREATE_INDEX=true

    # Remove duplicates
    #java -Xmx16g -jar /path/to/picard/picard.jar \
    #MarkDuplicates I=${ID}_mapped.bam \
    #O=rm_dup_${ID}_mapped.bam \
    #M=rm_dup_${ID}_mapped.metrics.txt \
    #REMOVE_DUPLICATES=true \
    #ASSUME_SORTED=true \
    #VALIDATION_STRINGENCY=LENIENT

    # Index the BAM file
    #samtools index rm_dup_${ID}_mapped.bam

    # Extract the non-duplicate reads
    #java -Xmx16g -jar /path/to/picard/picard.jar \
    #SamToFastq I=rm_dup_${ID}_mapped.bam \
    #F=rm_dup_${ID}_mapped_R1.fastq F2=rm_dup_${ID}_mapped_R2.fastq

    #gzip -f rm_dup_${ID}_mapped_R1.fastq
    #gzip -f rm_dup_${ID}_mapped_R2.fastq

    # Check XY non-par regions
    #samtools view -u \
    #../BAM/rm_dup_${ID}_mapped.bam | \
    #bedtools coverage -a nonpar.bed \
    #-b stdin -counts -F 0.1 | \
    #awk 'BEGIN{OFS="\t"}{print $0"\t"$4/($3-$2)"\t"$4}' | \
    #sed '1i CHR\tSTART\tEND\tBP\tDepth\tCount' > bedcov_${ID}_nonPAR_XY.txt
