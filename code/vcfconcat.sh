#!/bin/bash

# Description:
#   Concatenates VCF files split by chromosome. The input and output VCFs will have
# the same number of columns.
#
# Usage:
#  sbatch vcfconcat.sh

#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e

# Description:
#  Merge indexed vcf
# ignore chromose fo now. some merging problems by vcf-concat

source activate dropseq

module load vcftools

dir_pseudo="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489-pseudo"

cd $dir_pseudo

vcf-concat chr1.pseudo.sorted.vcf.gz chr10.pseudo.sorted.vcf.gz \
  chr11.pseudo.sorted.vcf.gz chr12.pseudo.sorted.vcf.gz \
  chr13.pseudo.sorted.vcf.gz chr14.pseudo.sorted.vcf.gz \
  chr15.pseudo.sorted.vcf.gz chr16.pseudo.sorted.vcf.gz \
  chr17.pseudo.sorted.vcf.gz chr18.pseudo.sorted.vcf.gz \
  chr19.pseudo.sorted.vcf.gz chr2.pseudo.sorted.vcf.gz \
  chr20.pseudo.sorted.vcf.gz chr21.pseudo.sorted.vcf.gz \
  chr22.pseudo.sorted.vcf.gz |\
  bgzip -c > pseudo.sorted.vcf.gz



#tabix -f -p vcf pseudo.sorted.vcf.gz


