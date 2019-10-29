#! /bin/bash

#
# SEVEN SCINTILLATING STEPS TO SPOTTING SIGNIFICANT SEQUENCE SUBSTITUTIONS
#
# (a.k.a., read to reference mapping with variant calling
#

# PLEASE NOTE: these commands are tailored for my specific data
# they will not work "out-of-the-box" for you!
# 
# you WILL have to change filenames
# you will almost certainly want to tweak arguments


# 1. index your reference genome
# this command creates a set of files that begin with the name of your reference genome file; e.g.:
# GCF_000018225.1_ASM1822v1_genomic.fna.amb
# GCF_000018225.1_ASM1822v1_genomic.fna.ann
# GCF_000018225.1_ASM1822v1_genomic.fna.bwt
# GCF_000018225.1_ASM1822v1_genomic.fna.pac
# GCF_000018225.1_ASM1822v1_genomic.fna.sa
# keep all of these files unchanged and together with your original reference genome (fna) file
bwa index GCF_000018225.1_ASM1822v1_genomic.fna


# 2. map your reads to your reference
# creates a sam file
bwa mem -t 12 -c 50000 -P -B 3 GCF_000018225.1_ASM1822v1_genomic.fna Driscol-M1-P1_S1_L001_R1_001.fastq.gz Driscol-M1-P1_S1_L001_R2_001.fastq.gz > S1.sam


# 3. convert your sam file to an ordered bam file
# an ordered bam file is required input for the next step (calling variants)
samtools view -b S1.sam | samtools sort -o S1.sorted.bam


# 4. index your sorted bam file
# keep the index file together with your sorted bam file
samtools index S1.sorted.bam S1.sorted.bai


# 5. call variants
# this creates a vcf file that contains information about every position in your reference genome
bcftools mpileup --threads 12 -6 -A -d 10000 -f GCF_000018225.1_ASM1822v1_genomic.fna S1.sorted.bam -O b -o S1.bcf
bcftools call -o S1.vcf -O v --threads 12 -A -m S1.bcf


# 6. filter your vcf file to extract just the rows that contain actual variants
egrep '\t\.\t[A-Z]+\t[A-Z]+\t' S1.vcf > S1.variants_only.vcf


# 7. combine your filtered vcf file with your reference genome (fna) and genomic feature (gff) files to identify the biological context of your variants
# this assume your input files are in the same dir as the vcf2table.pl script
./vcf2table.pl GCF_000018225.1_ASM1822v1_genomic.gff -g GCF_000018225.1_ASM1822v1_genomic.fna < S1.variants_only.vcf > S1.variant_table.txt


exit 1

