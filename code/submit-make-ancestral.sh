#!/bin/bash
set -e

chr_num=({2..22})
#chr_num=({1..3})
#chr_num=(1)

for i in ${chr_num[@]}; do

  sbatch make-ancestral.sh ${i}

done
