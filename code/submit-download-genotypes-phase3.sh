#!/bin/bash
set -e

http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

chr_num=({1..22} X)
#chr_num=(X Y)

for i in ${chr_num[@]}; do

  sbatch download-genotypes-phase3.sh ${website} ${i}

done
