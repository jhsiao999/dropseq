#!/bin/bash
set -e

#chr_num=({1..22} X Y)
chr_num=({1..22})
#chr_num=1

for i in ${chr_num[@]}; do

#  sbatch make-ancestral-scene-1.sh ${i}
  sbatch make-ancestral-scene-2.sh ${i}
  sbatch make-ancestral-scene-3.sh ${i}
  sbatch make-ancestral-scene-4.sh ${i}
#  sbatch make-ancestral-scene-5.sh ${i}

done
