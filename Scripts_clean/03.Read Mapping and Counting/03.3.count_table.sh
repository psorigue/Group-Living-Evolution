#!/bin/bash

# This script runs featureCounts for counting the reads mapped to the reference genome. 
# featureCounts version 2.0.2

# Paths and variables
PATH_FILES="${HOME}/03.Mapping_and_Counts/mapped_files/"
FILE_OUT="${HOME}/03.Mapping_and_Counts/counts_all_samples.txt"
FILE_ANNOT="${HOME}/Ref_genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_filtered_for_counts.gtf"
SAMPLE_ARRAY="${HOME}/sample_array.txt"
THR=8

# Read the sample names from text file
LIST_SAMPLES=( $( cat "${SAMPLE_ARRAY}" ) )

# Change directory to the path of the files
cd "${PATH_FILES}"

# Run featureCounts for all the samples together
featureCounts -O -T "${THR}" -a "${FILE_ANNOT}" -s 2 -o "${FILE_OUT}" $( for FILE in "${LIST_SAMPLES[@]}" ; do echo "${PATH_FILES}/${FILE}/${FILE}_Aligned.sortedByCoord.out.bam" ; done )
# The output dataset needs to be manually curated for downstream analyses.

exit