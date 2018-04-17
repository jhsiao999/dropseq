#!/bin/bash
set -e

source activate dropseq

module load vcftools

chr_num=({1..22})
#chr_num=(1 )

for i in ${chr_num[@]}; do

  sbatch get-genotypes-multiple.sh ${i}

done
