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

  rm ${DIR_OUT}/GATK/${ID}/Unfiltered_GATK_${ID}.vcf

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
  -O ${DIR_OUT}/GATK/${ID}/GATK_${ID}.filtered.vcf

  # Filter additionally for biallelic variants
  singularity exec ./bcftools.sif bcftools view -m2 -M2 -v snps -i "FILTER='PASS'" ${DIR_OUT}/GATK/${ID}/GATK_${ID}.filtered.vcf > ${DIR_OUT}/GATK/${ID}/GATK_${ID}.vcf
  rm ${DIR_OUT}/GATK/${ID}/GATK_${ID}.filtered.vcf

done
