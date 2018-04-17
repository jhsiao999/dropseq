#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e
website=$1
chr_num=$2

# download genotype data
url_genotypes=${website}/ALL.chr${chr_num}GRCh38.genotypes.20170504.vcf.gz

cd /project2/gilad/joycehsiao/dropseq-pipeline/genotypes

wget ${url_genotypes}

