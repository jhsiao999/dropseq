#!/bin/bash

# Description:
#   Concatenates VCF files split by chromosome. The input and output VCFs will have 
# the same number of columns.
#
# Usage:
#  sbatch vcf-concat.sh

#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e

# Description:
#  Merge indexed vcf
# ignore chromose fo now. some merging problems by vcf-concat

source activate dropseq

module load vcftools

dir_multiple="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-multiple"

cd $dir_multiple

vcf-concat chr1.genotypes.multiple.sorted.vcf.gz chr10.genotypes.multiple.sorted.vcf.gz \
  chr11.genotypes.multiple.sorted.vcf.gz chr12.genotypes.multiple.sorted.vcf.gz \
  chr13.genotypes.multiple.sorted.vcf.gz chr14.genotypes.multiple.sorted.vcf.gz \
  chr15.genotypes.multiple.sorted.vcf.gz chr16.genotypes.multiple.sorted.vcf.gz \
  chr17.genotypes.multiple.sorted.vcf.gz chr18.genotypes.multiple.sorted.vcf.gz \
  chr19.genotypes.multiple.sorted.vcf.gz chr2.genotypes.multiple.sorted.vcf.gz \
  chr20.genotypes.multiple.sorted.vcf.gz chr21.genotypes.multiple.sorted.vcf.gz \
  chr22.genotypes.multiple.sorted.vcf.gz |\
  bgzip -c > genotypes.multiple.sorted.vcf.gz

tabix -f -p vcf genotypes.multiple.sorted.vcf.gz


