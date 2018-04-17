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
dir_temp_scene="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489-tmp-scene-1"
dir_pseudo="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489-pseudo-scene-1"

chr_num=$1

echo "------------------- Processing chr "${chr_num}" -------------------"

# Create a pseudo ancestral genotype -----------------------
# Description:
# 1. Keep the snps satisfying the following condition
#  a. presence of ancestral allele
#  b. ancestral allele is the derived allele
# 2. Ancestral genotype is modified to be 0|1: one reference, one alterante

#vcftools --gzvcf $dir_input/chr${chr_num}.genotypes.18489.vcf.gz \
#  --recode --out $dir_temp_scene/chr${chr_num}.pseudo

# Extract snp positions for ancestral alelle genotyping

# 1. keep sites if ancestral allele present 
#cat $dir_temp/chr${chr_num}.genotypes.txt |
#  awk 'NR>1 {print}' |\
#  awk '{ if (!(($7!=".") OR ($7!="-"))) {print}}' \
#  > $dir_temp_scene/chr${chr_num}.genotypes.ancestral.1.txt

cat $dir_temp/chr${chr_num}.genotypes.txt |
  awk 'NR>1 {print}' |\
  awk '{ if (($7!=".") && ($7!="-"))  {print}}' \
  > $dir_temp_scene/chr${chr_num}.genotypes.ancestral.1.txt

# 2. keep sites if genotyped allele pair is not the same
# as the reference/alternate pair
cat $dir_temp_scene/chr${chr_num}.genotypes.ancestral.1.txt |\
  awk '{ if (!( ($3!=$5) && ($3!=$6) && ($4!=$5) && ($4!=$6) )) {print}}' |\
  awk 'BEGIN{print "#CHROM","POS","ALLELE1","ALLELE2","REF","ALT","AA"}1' \
  > $dir_temp_scene/chr${chr_num}.genotypes.ancestral.2.txt

# Extract position list
cat $dir_temp_scene/chr${chr_num}.genotypes.ancestral.2.txt | \
  awk '{print $1,$2}' > $dir_temp_scene/chr${chr_num}.genotypes.ancestral.positions.txt

# subset vcf
vcftools --gzvcf $dir_input/chr${chr_num}.genotypes.18489.vcf.gz \
  --positions $dir_temp_scene/chr${chr_num}.genotypes.ancestral.positions.txt \
  --recode --out $dir_temp_scene/chr${chr_num}.pseudo

# Make pseudo ancestral genotype and append to the vcf
cat $dir_temp_scene/chr${chr_num}.pseudo.recode.vcf | \
  awk 'NR>130 {print}' > \
  $dir_temp_scene/chr${chr_num}.pseudo.recode.main.vcf
cat $dir_temp_scene/chr${chr_num}.pseudo.recode.vcf | \
  awk 'NR<131 {print}' > \
  $dir_temp_scene/chr${chr_num}.pseudo.recode.header.vcf

# make pseudo ancestral genotype
# all is 0|1: reference and derived
cat $dir_temp_scene/chr${chr_num}.pseudo.recode.main.vcf | \
  awk 'NR>0 {print $11, "0|1"}' | \
  awk 'NR==1 {gsub("0\\|1","AA", $1);quit};1' \
  > $dir_temp_scene/chr${chr_num}.genotypes.ancestral.tmp.txt

# combine info
paste $dir_temp_scene/chr${chr_num}.pseudo.recode.main.vcf \
  $dir_temp_scene/chr${chr_num}.genotypes.ancestral.tmp.txt \
  > tmp && mv tmp $dir_temp_scene/chr${chr_num}.pseudo.tmp.vcf

cat $dir_temp_scene/chr${chr_num}.pseudo.recode.header.vcf \
  $dir_temp_scene/chr${chr_num}.pseudo.tmp.vcf |\
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








