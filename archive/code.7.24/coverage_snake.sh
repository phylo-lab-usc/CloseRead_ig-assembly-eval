#!/bin/bash
#SBATCH --job-name=coverage    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10                    # Run on a single CPU 
#SBATCH --time=6:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/coverage%j.log   # Standard output and error log
#SBATCH --mem=100G


source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly



while getopts s:a:b:f: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        a) assemblies=${OPTARG};;
        b) bam=${OPTARG};;
        f) loci=${OPTARG};;
    esac
done
HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality

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


touch "/home1/zhuyixin/zhuyixin_proj/AssmQuality/errorStats/${species}/pileup.end"

