#!/bin/bash

# Description
#  Parse sequence ID to get cell barcodes and molecules. 
#  Add cell barcodes and molecules to the bam file.

#SBATCH --time=36:00:00
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=6GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu

set -e

source activate dropseq

dir_bam="/project2/gilad/spott"
dir_tmp="/project2/gilad/joycehsiao/dropseq-pipeline/genotypes-human-control-tmp"
dir_output="/project2/gilad/joycehsiao/dropseq-pipeline/bam-indexed-human-control"

# Extract cell barcodes
samtools view $dir_bam/Aligned.SortedByCoordinate.out.bam |\
  awk '{split($1, a, "_"); print a[2]}' > $dir_tmp/barcodes.txt

cat $dir_tmp/barcodes.txt |\
  awk '{if (NR > 0) {print "CB:Z:"$0}}' > tmp && mv tmp $dir_tmp/barcodes.txt

# Extract molecules
samtools view $dir_bam/Aligned.SortedByCoordinate.out.bam |\
  awk '{split($1, a, "_"); print a[3]}' > $dir_tmp/molecules.txt
cat $dir_tmp/molecules.txt |\
  awk '{if (NR > 0) {print "UB:Z:"$0}}' > tmp && mv tmp $dir_tmp/molecules.txt

# Export bam header
samtools view -H $dir_bam/Aligned.SortedByCoordinate.out.bam \
  > $dir_tmp/header.sam

# Convert bam to sam
samtools view $dir_bam/Aligned.SortedByCoordinate.out.bam > $dir_tmp/Aligned.SortedByCoordinate.out.sam

# Add barcode and molecule field to sam
paste $dir_tmp/Aligned.SortedByCoordinate.out.sam $dir_tmp/barcodes.txt $dir_tmp/molecules.txt \
  > $dir_tmp/tmp.sam && mv $dir_tmp/tmp.sam $dir_tmp/Aligned.SortedByCoordinate.out.sam

# Add hearder field to sam
cat $dir_tmp/header.sam $dir_tmp/Aligned.SortedByCoordinate.out.sam  \
  > $dir_tmp/tmp && mv $dir_tmp/tmp $dir_output/Aligned.SortedByCoordinate.out.sam

# Convert sam to bam
samtools view -bS $dir_output/Aligned.SortedByCoordinate.out.sam |\
  samtools sort -O bam -o $dir_output/Aligned.SortedByCoordinate.out.bam

#samtools index tmp3_sorted.bam tmp3_sorted.bai


