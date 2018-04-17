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
dir_output="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-18489"

vcftools --gzvcf $dir_input/ALL.chr${chr_num}_GRCh38.genotypes.20170504.vcf.gz \
  --indv NA18489 \
  --remove-indels \
  --recode --recode-INFO-all \
  --stdout | \
  gzip -c > $dir_output/chr${chr_num}.genotypes.18489.vcf.gz



