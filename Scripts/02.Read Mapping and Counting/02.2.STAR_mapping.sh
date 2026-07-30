#!/bin/bash

# This script runs STAR for mapping the reads to the reference genome.
# STAR version 2.7.11b
# Samtools version 1.21

# Paths and variables
PATH_FILES="${HOME}/raw_data/fastp/" # Path to the fastq files after running fastp
PATH_OUT="${HOME}/03.Mapping_and_Counts/mapped_files/"
PATH_GENOME="${HOME}/03.Mapping_and_Counts/STAR_genome/"
REFERENCE_GTF="${HOME}/Ref_genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf"
SAMPLE_ARRAY="${HOME}/sample_array.txt"
THR=8

# Read the sample names from text file
LIST_SAMPLES=( $( cat "${SAMPLE_ARRAY}" ) )

# Change directory to the path of the files
cd "${PATH_FILES}"

# Loop through the samples and run STAR for each sample
for FILE in "${LIST_SAMPLES[@]}" ; do

    # Create output directory for each sample
    mkdir -p "${PATH_OUT}/${FILE}"
    cd "${PATH_OUT}/${FILE}"
    
    # Run STAR for each sample. 
    STAR --genomeDir "${PATH_GENOME}" \
        --sjdbGTFfile "${REFERENCE_GTF}" \
        --readFilesIn "${PATH_FILES}"/"${FILE}".out.fastq \
        --outFileNamePrefix "./${FILE}_" \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN "${THR}" \
        --outFilterMultimapNmax 10 \
        --outFilterMatchNminOverLread 0.4 \
        --outFilterScoreMinOverLread 0.4 \
        --quantMode TranscriptomeSAM \
        --alignIntronMax 100000
    
    # Index the BAM file
    samtools index "./${FILE}_*.bam"
    
    cd "${PATH_OUT}"
    
done

exit