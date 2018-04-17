#!/bin/bash

# Description:
#   Make ancestral genotype VCF
# for one chromosome
#
# Usage:
#   bash preprocess-ancestral.sh ${chr_num}

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

chr_num=$1

echo "------------------- Processing chr "${chr_num}" -------------------"

# Extract AA info from the genotype file -----------------------------------
# print the AA field, replace header, remove "AA=." from the content
zcat $dir_input/chr${chr_num}.genotypes.18489.vcf.gz |\
  awk 'NR>131 {print $8}'   | \
    awk -F";" '{for(i=1;i<=NF;i++){if ($i ~ /AA/){print $i}}}' | \
    awk 'NR>0 {gsub("AA=", "", $0); gsub("\\|", "", $0); quit};1' \
    > $dir_temp/chr${chr_num}.AA.txt

# Extract reference/alternate info from the genotype file ----------------
# print the AA field, replace header, remove "AA=." from the content
zcat $dir_input/chr${chr_num}.genotypes.18489.vcf.gz | \
  awk 'NR>131 {print $4, $5}' > $dir_temp/chr${chr_num}.ref.alt.txt

# Get allele counts --------------------------------------
# Use vcftools to output allele counts
vcftools --gzvcf \
  $dir_input/chr${chr_num}.genotypes.18489.vcf.gz \
  --counts --out $dir_temp/chr${chr_num}

# Edit allele count file
cat $dir_temp/chr${chr_num}.frq.count | \
  awk 'NR>1 {print $1,$2,$5,$6}' | \
# process the first allele
  awk '{gsub("A:0", "N", $3); gsub("C:0", "N", $3); gsub("G:0", "N", $3); gsub("T:0", "N", $3); quit};1' | \
  awk '{gsub("A:1", "A", $3); gsub("C:1", "C", $3); gsub("G:1", "G", $3); gsub("T:1", "T", $3); quit};1' | \
  awk '{gsub("A:2", "A", $3); gsub("C:2", "C", $3); gsub("G:2", "G", $3); gsub("T:2", "T", $3); quit};1' | \
# process the second allele
  awk '{gsub("A:0", "N", $4); gsub("C:0", "N", $4); gsub("G:0", "N", $4); gsub("T:0", "N", $4); quit};1' | \
  awk '{gsub("A:1", "A", $4); gsub("C:1", "C", $4); gsub("G:1", "G", $4); gsub("T:1", "T", $4); quit};1' | \
  awk '{gsub("A:2", "A", $4); gsub("C:2", "C", $4); gsub("G:2", "G", $4); gsub("T:2", "T", $4); quit};1'  \
  > $dir_temp/chr${chr_num}.alleles.txt

# Combine genotype information -----------------------------
paste $dir_temp/chr${chr_num}.alleles.txt $dir_temp/chr${chr_num}.ref.alt.txt |\
  column -s $'\t' -t > $dir_temp/chr${chr_num}.genotypes.tmp.txt
paste $dir_temp/chr${chr_num}.genotypes.tmp.txt $dir_temp/chr${chr_num}.AA.txt | \
  awk 'BEGIN{print "#CHROM","POS","ALLELE2","ALLELE2","REF","ALT","AA"}1' \
  > tmp.${chr_num} && mv tmp.${chr_num} $dir_temp/chr${chr_num}.genotypes.txt



