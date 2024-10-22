#!/bin/sh
<<VARIABLE
REF_FASTA: reference genome file (ex: hg37_1kg_decoy)
DBSNP_FILE: dbsnp file
VARIABLE

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

# Step 1: Variant calling and filtering for all samples
for ID in ${SAMPLES}; do
  mkdir -p ${DIR_OUT}/GATK/${ID}
  cd ${DIR_OUT}/GATK/${ID}

  # Call variants with GATK HaplotypeCaller
  singularity exec gatk-4.sif gatk \
  --java-options "-Xmx12g" \
  HaplotypeCaller \
  -R ${REF_FASTA} \
  -I ${DIR_IN}/${ID}_rg_nodup.bam \
  --dbsnp ${DBSNP_FILE} \
  -O Unfiltered_GATK_${ID}.vcf

  # Filter to retain SNPs and autosomes
  for i in `seq 1 22`; do echo -e ${i}; done > auto.list

  singularity exec gatk-4.sif gatk \
  --java-options "-Xmx12g" \
  SelectVariants \
  -V Unfiltered_GATK_${ID}.vcf \
  -select-type SNP \
  -L auto.list \
  -O Auto_Unfiltered_GATK_${ID}.vcf

  rm auto.list
  rm Unfiltered_GATK_${ID}.vcf

  # Apply hard filters to SNP dataset
  singularity exec gatk-4.sif gatk \
  --java-options "-Xmx8g" \
  VariantFiltration \
  -V Auto_Unfiltered_GATK_${ID}.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "DP < 2.0" --filter-name "DP2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O GATK_${ID}.filtered.vcf

  singularity exec ./bcftools.sif bcftools view -i "%FILTER='PASS'" GATK_${ID}.filtered.vcf > GATK_${ID}.vcf
  rm GATK_${ID}.filtered.vcf

  # Step 2: Convert each filtered VCF to plink format
  singularity exec plink.sif plink1.9 \
  --vcf GATK_${ID}.vcf \
  --make-bed \
  --out ${ID}_plink
done

# Step 3: Merge all sample plink datasets into one
# Create a file containing the plink dataset file names
plink_files=""
for ID in ${SAMPLES}; do
  plink_files="${plink_files} ${ID}_plink"
done

# If more than one sample, merge them
if [ $(echo ${SAMPLES} | wc -w) -gt 1 ]; then
  echo ${plink_files} | tr ' ' '\n' > merge_list.txt
  plink --bfile $(echo ${plink_files} | awk '{print $1}') --merge-list merge_list.txt --make-bed --out merged_samples
else
  cp ${ID}_plink.bed merged_samples.bed
  cp ${ID}_plink.bim merged_samples.bim
  cp ${ID}_plink.fam merged_samples.fam
fi

# Step 4: Perform multi-sample genotype comparison with plink
# Calculate identity-by-state (IBS) or identity-by-descent (IBD)
singularity exec plink.sif plink1.9 \
--bfile merged_samples \
--genome \
--out sample_comparison

# The --genome command will output a pairwise comparison of IBS/IBD between all samples.
# Use this output to determine which samples are genetically similar and visualize the similarity matrix
