#!/bin/bash

# Description:
#   Make ancestral genotype VCF
# for one chromosome
#
# Usage:
#   bash make-ancestral.sh ${chr_num}

#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e

source activate dropseq

dir_input="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489"
dir_temp="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489-tmp"
dir_pseudo="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489-pseudo"

chr_num=$1

# Create a pseudo ancestral genotype -----------------------
# Description:
# 1. Keep the snps satisfying the following condition
#  a. presence of ancestral allele
#  b. ancestral allele is the derived allele
# 2. Ancestral genotype is modified to be 0|1: one reference, one alterante


echo "------------------- Processing chr "${chr_num}" -------------------"

# Output SNP positions, reference/alternate alleles, and individual genotype information
vcftools --gzvcf $dir_input/chr${chr_num}.genotypes.18489.vcf.gz \
  --recode --out $dir_temp/chr${chr_num}.pseudo

# keep sites if ancestral allele present
cat $dir_temp/chr${chr_num}.genotypes.txt |
  awk 'NR>1 {print}' |\
  awk '{ if (($7!=".") OR ($7!="-")) {print}}' \
  > $dir_temp/chr${chr_num}.genotypes.ancestral.1.txt

# approach 1: keep those that ancestral allele is either reference or alternate allele
# then assign 0|0 for ancestral alleles matching reference alleles,
# and 1|1 for ancestral alleles matching alternate alleles

# ancestral matches reference (0|0)
cat $dir_temp/chr${chr_num}.genotypes.ancestral.1.txt |\
  awk '{print toupper($0)}' | \
  awk '{ if ($5==$7) {print}}' | \
  awk '{print $1,$2}' | \
  awk 'NR>0 {print $0, "0|0"}' > $dir_temp/chr${chr_num}.genotypes.ancestral.2.txt

# ancestral matching altenate (1|1)
cat $dir_temp/chr${chr_num}.genotypes.ancestral.1.txt |\
  awk '{print toupper($0)}' | \
  awk '{ if ($6==$7) {print}}' | \
  awk '{print $1,$2}' | \
  awk 'NR>0 {print $0, "1|1"}' > $dir_temp/chr${chr_num}.genotypes.ancestral.3.txt

# Combine SNPs
cat $dir_temp/chr${chr_num}.genotypes.ancestral.2.txt \
  $dir_temp/chr${chr_num}.genotypes.ancestral.3.txt > \
  $dir_temp/chr${chr_num}.genotypes.ancestral.positions.txt

sort -g -k2 $dir_temp/chr${chr_num}.genotypes.ancestral.positions.txt > \
  tmp && mv tmp $dir_temp/chr${chr_num}.genotypes.ancestral.positions.txt

cat $dir_temp/chr${chr_num}.genotypes.ancestral.positions.txt | \
  awk '{print $1,$2}' > $dir_temp/chr${chr_num}.genotypes.ancestral.positions.list.txt

# Subset vcf to include snps with ancestral allele present
vcftools --gzvcf $dir_input/chr${chr_num}.genotypes.18489.vcf.gz \
  --positions $dir_temp/chr${chr_num}.genotypes.ancestral.positions.list.txt \
  --recode --out $dir_temp/chr${chr_num}.pseudo

# Make pseudo ancestral genotype (homozygous) and append to the vcf
# first separate header from the vcf body
cat $dir_temp/chr${chr_num}.pseudo.recode.vcf | \
  awk 'NR>130 {print}' > \
  $dir_temp/chr${chr_num}.pseudo.recode.main.vcf
cat $dir_temp/chr${chr_num}.pseudo.recode.vcf | \
  awk 'NR<131 {print}' > \
  $dir_temp/chr${chr_num}.pseudo.recode.header.vcf

# make pseudo ancestral genotype
cat $dir_temp/chr${chr_num}.genotypes.ancestral.positions.txt  | \
  awk 'NR>0 {print $3}' | \
  awk 'BEGIN{print "AA"}; {print};' \
  > $dir_temp/chr${chr_num}.genotypes.ancestral.tmp.txt

# combine info
paste $dir_temp/chr${chr_num}.pseudo.recode.main.vcf \
  $dir_temp/chr${chr_num}.genotypes.ancestral.tmp.txt \
  > tmp && mv tmp $dir_temp/chr${chr_num}.pseudo.tmp.vcf

cat $dir_temp/chr${chr_num}.pseudo.recode.header.vcf \
  $dir_temp/chr${chr_num}.pseudo.tmp.vcf |\
  bgzip -c > $dir_pseudo/chr${chr_num}.pseudo.vcf.gz

# Prepare VCF for mapping ---------------------
# bam is ucsc chromosome name; vcf is Ensembl chromosome name
# add "chr" to chromosome name
zcat $dir_pseudo/chr${chr_num}.pseudo.vcf.gz |\
  awk '{if (NR > 131) {print "chr"$0}}' |\
  gzip -c > $dir_pseudo/chr${chr_num}.pseudo.main.vcf.gz

zcat $dir_pseudo/chr${chr_num}.pseudo.vcf.gz |\
  awk '{if (NR <= 131) {print}}' |\
  gzip -c > $dir_pseudo/chr${chr_num}.pseudo.header.vcf.gz

gunzip -f $dir_pseudo/chr${chr_num}.pseudo.main.vcf.gz
gunzip -f $dir_pseudo/chr${chr_num}.pseudo.header.vcf.gz

cat $dir_pseudo/chr${chr_num}.pseudo.header.vcf $dir_pseudo/chr${chr_num}.pseudo.main.vcf|\
  bgzip -c -f > $dir_pseudo/chr${chr_num}.pseudo.sorted.vcf.gz

tabix -f -p vcf $dir_pseudo/chr${chr_num}.pseudo.sorted.vcf.gz

rm $dir_pseudo/chr${chr_num}.pseudo.main.vcf
rm $dir_pseudo/chr${chr_num}.pseudo.header.vcf








