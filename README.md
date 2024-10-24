# Analysis of the human reads in gut metagenome shotgun sequencing data
This is a fork of the repository of the codes used in the Tomofuji et al (Reconstruction of the personal information from the contaminated human reads in the gut metagenome shotgun sequencing data).  
Code is simplified to recover only following information from the human reads in the metagenome shotgun sequencing data: <br>
・Genetic sex of the metagenome shotgun sequencing data <br>
・Call variants so that SNPs can be used to identify samples from the same person in a dataset (relaxed criteria for low coverage) <br>

# Prepare for running pipelines 1a, 1b, 5a

## Software
・singularity, build .sif files:

```
# One container for cleaning reads, extracting human reads (masked overlap with microbial similar regions)
singularity build cleanreads.sif docker://chvav/cleanreads@sha256:b8803f1bdb8a91b76286a5d26b38c9b068120869f6bf829c000d3364a6933359

# Bowtie 2 v2.3.4.1 and samtools v1.7
singularity build bowtie2_samtools.sif docker://davelabhub/bowtie2_plus_samtools@sha256:04afb9762780dcac25ebd1947c2e12ef81f96a51604dcaed5fb2081804051108

# Picard v2.21.5
singularity build ./picard.sif docker://pegi3s/picard@sha256:927475a7ade08cb9376127f04d410a9c78bf7f4fd9a6791363ceed68cd8aac74

# bedtools v2.30.0
singularity build ./bedtools.sif docker://pegi3s/bedtools@sha256:fc96153f83e4c69d742ee9b2bde33d5d297ef83b0a9b2d81ae32e14e26deb8fc

# htslib 1.21
singularity build htslib.sif docker://staphb/htslib@sha256:b6af37ebf92852d4cace4d9e9f1dd434e93245800f046b7788695fcdc8418b8c

# GATK 4.1.4.1 (Picard included here as well)
singularity build ./gatk-4.sif docker://pegi3s/gatk-4@sha256:8415200fe87b4d4eb295e9decb207727e9f75aad73c860a8d0ffe3385a8eaa7a

# bcftools 1.21
sudo singularity build ./bcftools.sif docker://staphb/bcftools@sha256:c256fe5b7f8a216c517c2c6ff8fd68a6643ee614e93eedd181aa2c53600ec209

```

## Files

Download the bbduk provided file with Illumina adaptor and contaminant sequences (PhiX spike in control):
```
wget https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa
```

Download the version human reference genome, masked for regions that may be similar to microbial genomes (https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/37175-introducing-removehuman-human-contaminant-removal#post245912):
https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?resourcekey=0-PsIKmg2q4EvTGWGOUjsKGQ

This file will be used for initial extraction of the human reads.

Index the file:
```
singularity exec cleanreads.sif bbmap.sh ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
```

Download the human reference genome and index for final mapping:

```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
singularity exec bowtie2_samtools.sif bowtie2-build human_g1k_v37.fasta.gz human_g1k_v37
```
Save a nonpar.bed file (3 columns, tab-seperated) with the non-pseudoautosomal regions of the XY chromosomes (corresponding to GRCh37/hg19): <br>
X       2699521 154931043 <br>
Y       2649521 59034049 <br>

Download the corresponding SNP database from GATK resource bundle:
```
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
```

Prepare bgzip file versions of the human genome and corresponding snp database, create indices and sequence dictionary:

```
# Prep human genome
gunzip -k human_g1k_v37.fasta.gz
singularity exec htslib.sif bgzip -c human_g1k_v37.fasta > human_g1k_v37_bgzip.fasta.gz
singularity exec bowtie2_samtools.sif samtools faidx human_g1k_v37_bgzip.fasta.gz

singularity exec gatk-4.sif gatk CreateSequenceDictionary \
  -R human_g1k_v37_bgzip.fasta.gz \
  -O human_g1k_v37_bgzip.dict

# Prep snp db
gunzip -k dbsnp_138.b37.vcf.gz
singularity exec htslib.sif bgzip -c dbsnp_138.b37.vcf > dbsnp_138.b37_bgzip.gz
singularity exec gatk-4.sif gatk IndexFeatureFile -I dbsnp_138.b37_bgzip.gz
```

# 1. Extraction of the human reads and prediction of genetic xex

## Clean reads and extract human reads

```
chmod u+x PIPELINE_1a_cleanreads.sh
./PIPELINE_1a_cleanreads.sh <DIR_IN> <DIR_OUT>
```

## Predict sex

```
chmod u+x PIPELINE_1b_predict_sex.sh
./PIPELINE_1b_predict_sex.sh <DIR_IN> <DIR_OUT>
```

## Call SNPs

```
chmod u+x PIPELINE_5a_variant_call_GATK.sh
./PIPELINE_5a_variant_call_GATK.sh <DIR_IN> <DIR_OUT>
```
