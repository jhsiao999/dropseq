#!/bin/bash
set -e

source activate dropseq

module load vcftools

chr_num=({1..22})
#chr_num=(1 2)

for i in ${chr_num[@]}; do

  sbatch get-genotypes-18489.sh ${i}

done
