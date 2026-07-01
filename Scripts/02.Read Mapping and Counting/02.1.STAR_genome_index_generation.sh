#!/bin/bash

# This script generates the STAR genome index for the reference genome. This step is necessary before running STAR for mapping the reads to the reference genome.
# STAR version 2.7.11b

# Paths and variables
PATH_GENOME_GENERATION="${HOME}/03.Mapping_and_Counts/STAR_genome/"
REF_FASTA="${HOME}/Ref_genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa"
REF_GTF="${HOME}/Ref_genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf"
THR=8

# Change directory to the path of the genome generation
cd "${PATH_GENOME_GENERATION}"

# Generate the STAR genome index
STAR --runThreadN "${THR}" --runMode genomeGenerate --genomeDir "${PATH_GENOME_GENERATION}" --genomeFastaFiles "${REF_FASTA}" --sjdbGTFfile "${REF_GTF}" --genomeSAindexNbases 13 --sjdbOverhang 99

exit
