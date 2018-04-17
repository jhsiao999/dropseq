#!/bin/bash

# 1. Download genotypes and snps for hg38(GRCh38) in 1000 Genomes Project.
#
# Input:
#   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
#
# Output:
#   genotypes/
#
# Description:
#    Note that these are based on dbSNP149 and a "liftover" mapping from the GRCh37
# genome assembly used by the 1000 Genomes Project to the newer GRCh38 assembly.
# Link to 1000 Genomes Project: http://www.internationalgenome.org/announcements/
#
# submit-download-genotypes downloads doata using download-genotypes.sh
bash submit-download-genotypes.sh


# 2. Extract genotypes for the individual included in the experiment
#    from 1000 Genome Project Phase 3
#
# Input:
#   genotypes/
# Output:
#   genotypes-18489/
bash get-genotypes-18489.sh


# 3. Preprocess bam files
#
# Description:
#   Modify bam files to add UMI and cell barcode fields,
# because demuxlet requires UMI and cell barcode fields.
#
# a. for data from human-chip mix
#
# Input:
#   /project2/gilad/spott/Dropseq/iPSC_HS_panTro/data/sorted_reads
# Output:
#   bam_indexed/
bash parse-bam-index.sh

# b. for data from a single human sample
#
# Input:
#   /project2/gilad/spott
#
# Output:
#   bam-indexed-human-control/
bash parse-bam-index-human-control.sh



# 4. Prepare allele information
#
# Description:
#
#   Make a table that contains individual genotypes,
# reference/alternative genotypes, and ancestral genotypes.
# This table is used as reference for subsetting VCF files.
#
# Input:
#   genotypes-18489/
# Output:
#  genotypes-18489-tmp/
#  genotypes-18489-tmp/chr${chr_num}.genotypes.txt contains chromosome-level data
#
# Usage:
#   submit-preprocess-ancestral submits preprocess-ancestral for one chromosome at a time
bash submit-preprocess-ancestral.sh

# This step has changed!
# much of the code is now run in make-ancestral-hetero.R
# 5. Make psuedo genotype for ancestral allele
#
# Input:
#   genotypes-18489/
# Output
#   genotypes-18489-pseudo/
# Usage:
#   submit-make-ancestral.sh submits make-ancestral.sh for one chromosome at a time
#bash submit-make-ancestral.sh
bash make-ancestral-hetero.sh

# 6. Combine VCF files
#
# Input:
#   genotypes-18489-pseudo/
# Output:
#   genotypes-18489-pseudo/
bash vcfconcat.sh

# 6. Apply demuxlet
# Input:
# Output
bash demuxlet.sh


