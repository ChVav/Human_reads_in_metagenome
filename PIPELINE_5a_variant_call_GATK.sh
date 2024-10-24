#!/bin/sh

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <DIR_IN> <DIR_OUT>"
    exit 1
fi

# Input and output directories from command line arguments
DIR_IN="$1"
DIR_OUT="$2"

# Create output directory if it doesn't exist
mkdir -p ${DIR_OUT}/GATK/PLINK_COMPARISON

# Extract sample IDs from BAM file names in DIR_IN
SAMPLES=$(find ${DIR_IN} -name '*_rg_nodup.bam' | sed -r 's!.*/(.*)_rg_nodup\.bam!\1!')

# Create interval list to retain SNPs on autosomes
for i in `seq 1 22`; do echo "${i}"; done > ${DIR_OUT}/GATK/${ID}/auto.list

# Step 1: Variant calling and filtering for all samples
for ID in ${SAMPLES}; do
  mkdir -p ${DIR_OUT}/GATK/${ID}
 
  # Call variants with GATK HaplotypeCaller
  singularity exec gatk-4.sif gatk \
  --java-options "-Xmx12g" \
  HaplotypeCaller \
  -R human_g1k_v37_bgzip.fasta.gz \
  -I ${DIR_IN}/${ID}_rg_nodup.bam \
  --dbsnp dbsnp_138.b37_bgzip.gz \
  -O ${DIR_OUT}/GATK/${ID}/Unfiltered_GATK_${ID}.vcf

  singularity exec gatk-4.sif gatk \
  --java-options "-Xmx12g" \
  SelectVariants \
  -V ${DIR_OUT}/GATK/${ID}/Unfiltered_GATK_${ID}.vcf \
  -select-type SNP \
  -L ${DIR_OUT}/GATK/auto.list \
  -O ${DIR_OUT}/GATK/${ID}/Auto_Unfiltered_GATK_${ID}.vcf

  #rm ${DIR_OUT}/GATK/${ID}/Unfiltered_GATK_${ID}.vcf

  # Apply hard filters to SNP dataset
  # MQRankSum and ReadPosRankSum not calculated (due to low coverage?)
  singularity exec gatk-4.sif gatk \
  --java-options "-Xmx8g" \
  VariantFiltration \
  -V ${DIR_OUT}/GATK/${ID}/Auto_Unfiltered_GATK_${ID}.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "DP < 2.0" --filter-name "DP2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  #-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  #-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O ${DIR_OUT}/GATK/${ID}/GATK_${ID}.filtered.vcf

  singularity exec ./bcftools.sif bcftools view -i "FILTER='PASS'" ${DIR_OUT}/GATK/${ID}/GATK_${ID}.filtered.vcf > ${DIR_OUT}/GATK/${ID}/GATK_${ID}.vcf
  #rm ${DIR_OUT}/GATK/${ID}/GATK_${ID}.filtered.vcf

  # Step 2: Convert each filtered VCF to plink format
  singularity exec plink.sif plink1.9 \
  --vcf ${DIR_OUT}/GATK/${ID}/GATK_${ID}.vcf \
  --make-bed \
  --out ${DIR_OUT}/GATK/PLINK_COMPARISON/${ID}_plink
done

# Step 3: Merge all sample plink datasets into one
# Create a file containing the plink dataset file names
plink_files=""
for ID in ${SAMPLES}; do
  plink_files="${plink_files} ${DIR_OUT}/GATK/PLINK_COMPARISON/${ID}_plink"
done

# If more than one sample, merge them
if [ $(echo ${SAMPLES} | wc -w) -gt 1 ]; then
  echo ${plink_files} | tr ' ' '\n' > ${DIR_OUT}/GATK/PLINK_COMPARISON/merge_list.txt
  singularity exec plink.sif plink1.9 --bfile $(echo ${plink_files} | awk '{print $1}') --merge-list ${DIR_OUT}/GATK/PLINK_COMPARISON/merge_list.txt --make-bed --out ${DIR_OUT}/GATK/PLINK_COMPARISON/merged_samples
else
  cp ${DIR_OUT}/GATK/PLINK_COMPARISON/${ID}_plink.bed ${DIR_OUT}/GATK/PLINK_COMPARISON/merged_samples.bed
  cp ${DIR_OUT}/GATK/PLINK_COMPARISON/${ID}_plink.bim ${DIR_OUT}/GATK/PLINK_COMPARISON/merged_samples.bim
  cp ${DIR_OUT}/GATK/PLINK_COMPARISON/${ID}_plink.fam ${DIR_OUT}/GATK/PLINK_COMPARISON/merged_samples.fam
fi

# Step 4: Perform multi-sample genotype comparison with plink
# Calculate identity-by-state (IBS) or identity-by-descent (IBD)
singularity exec plink.sif plink1.9 \
--bfile ${DIR_OUT}/GATK/PLINK_COMPARISON/merged_samples \
--genome \
--out ${DIR_OUT}/GATK/PLINK_COMPARISON/sample_comparison

# The --genome command will output a pairwise comparison of IBS/IBD between all samples.
# Use this output to determine which samples are genetically similar and visualize the similarity matrix
