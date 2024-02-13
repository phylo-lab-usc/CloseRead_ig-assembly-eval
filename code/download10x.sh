#!/bin/bash

# Base directory where subdirectories will be created for each sample
BASE_DIR="/home1/zhuyixin/sc1/AssmQuality/10x_fastq"

export BASE_DIR

download_file() {
    sample=$1
    url=$2
    directory="${BASE_DIR}/${sample}"
    mkdir -p "${directory}" # Ensure the directory exists
    filename=$(basename "${url}")
    # Use wget or curl to download the file
    wget -q -O "${directory}/${filename}" "${url}" || curl -s -o "${directory}/${filename}" "${url}"
    echo "Downloaded ${url} to ${directory}/${filename}"
}

export -f download_file

cat /home1/zhuyixin/sc1/AssmQuality/10x.url.txt | parallel -j 20 --colsep ' ' download_file {1} {2}

