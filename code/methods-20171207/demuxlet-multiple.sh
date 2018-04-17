#!/bin/bash

# Usage
#  sbatch demuxlet.sh

#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e

# Options
#  --sam : Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed.
#  --vcf: Input a single VCF/BCF file, containing the individual genotypes (GT), 
#         posterior probability (GP), or genotype likelihood (PL).
#  --field: FORMAT field to extract the genotype, likelihood, or posterior from.
#  --out: Output file prefix

source activate dropseq

# path to demuxlet
demuxlet="/project/gilad/software/midway2/demuxlet/bin/demuxlet"

dir_bam_indexed="/project2/gilad/joycehsiao/dropseq-pipeline/bam_indexed"
dir_pseudo="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-multiple"
dir_output="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-multiple-output"

$demuxlet --sam $dir_bam_indexed/assigned_sorted.bam \
  --vcf $dir_pseudo/genotypes.multiple.sorted.vcf.gz \
  --field GT --out $dir_output/genotypes.multiple

