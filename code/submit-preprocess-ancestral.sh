#!/bin/bash
set -e

#chr_num=({1..22} X Y)
chr_num=({1..22})
#chr_num=(1)

for i in ${chr_num[@]}; do

  sbatch preprocess-ancestral.sh ${i}

done
