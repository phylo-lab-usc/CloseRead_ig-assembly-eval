#!/bin/bash

# Specify the directory where you want to save the FASTQ files
SRA_DIR="/home1/zhuyixin/zhuyixin_proj/AssmQuality/hifi_fastq"

# Function to convert a single SRA file to FASTQ
convert_sra_to_fastq() {
    SRA_FILE=$1
    echo $SRA_FILE >&2
    fastq-dump  --split-files --gzip "${SRA_FILE}"
    echo "Converted ${SRA_FILE} to FASTQ" >&2
}

export -f convert_sra_to_fastq

# Find all .sra files and pass them to parallel for conversion
find "${SRA_DIR}" -name "*.sra" | parallel -j 4 convert_sra_to_fastq

echo "All downloads completed."
