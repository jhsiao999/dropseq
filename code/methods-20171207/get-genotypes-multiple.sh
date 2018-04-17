#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e
chr_num=$1

echo "process chr"${chr_num}

dir_input="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes"
dir_output="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-multiple"
dir_tmp="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-multiple-tmp"

vcftools --gzvcf $dir_input/ALL.chr${chr_num}_GRCh38.genotypes.20170504.vcf.gz \
  --indv NA18489 --indv NA18499 --indv NA18501 --indv NA18502 --indv NA18507 --indv NA18508\
  --remove-indels \
  --recode --recode-INFO-all \
  --stdout | \
  bgzip -c > $dir_output/chr${chr_num}.genotypes.multiple.vcf.gz

#tabix -f -p vcf $dir_output/chr${chr_num}.genotypes.multiple.vcf.gz


# bam is ucsc chromosome name; vcf is Ensembl chromosome name
# add "chr" to chromosome name
zcat $dir_output/chr${chr_num}.genotypes.multiple.vcf.gz |\
  awk '{if (NR > 131) {print "chr"$0}}' |\
  gzip -c > $dir_tmp/chr${chr_num}.main.vcf.gz

zcat $dir_output/chr${chr_num}.genotypes.multiple.vcf.gz |\
  awk '{if (NR <= 131) {print}}' |\
  gzip -c > $dir_tmp/chr${chr_num}.header.vcf.gz

gunzip -f $dir_tmp/chr${chr_num}.main.vcf.gz
gunzip -f $dir_tmp/chr${chr_num}.header.vcf.gz

cat $dir_tmp/chr${chr_num}.header.vcf $dir_tmp/chr${chr_num}.main.vcf |\
  bgzip -c -f > $dir_output/chr${chr_num}.genotypes.multiple.sorted.vcf.gz

tabix -f -p vcf $dir_output/chr${chr_num}.genotypes.multiple.sorted.vcf.gz

rm $dir_tmp/chr${chr_num}.main.vcf
rm $dir_tmp/chr${chr_num}.header.vcf

