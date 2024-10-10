#!/bin/sh
while getopts s:a:b:f:d:c:e: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        a) assemblies=${OPTARG};;
        b) bam=${OPTARG};;
        f) loci=${OPTARG};;
        d) HOME=${OPTARG};;
        c) conda=${OPTARG};;
        e) condaEnv=${OPTARG};;
    esac
done

source ${conda}
conda init
conda activate ${condaEnv}/ig-assembly-eval

echo ${bam}

# Loop through each line of the file
while IFS= read -r line; do
    # Split the line into fields
    read -ra arr <<< "$line"
    # Extract information
    species=${arr[0]}
    haplotype=${arr[1]}
    gene=${arr[2]}
    loc=${arr[3]}:${arr[4]}-${arr[5]}

    # Check if the haplotype is primary
    if [ "$haplotype" = "primary" ]; then
        # Construct the output file path
        output_file="${HOME}/errorStats/${species}/${gene}_pri_pileup.txt"
        # Ensure output directory exists
        mkdir -p "$(dirname "$output_file")"
    else
        output_file="${HOME}/errorStats/${species}/${gene}_alt_pileup.txt"
        # Ensure output directory exists
        mkdir -p "$(dirname "$output_file")"
    fi
    samtools mpileup -Q 0 -q 0 -aa -f "${assemblies}" -r "$loc" "${bam}" >> "$output_file"
done < "$loci"


touch "${HOME}/errorStats/${species}/pileup.end"

