#!/bin/bash
#SBATCH --job-name=cigar_processing
#SBATCH --array=1-4%4
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --output=log/cigar_processing_%A_%a.log

# Define the thresholds
THRESHOLDS=(0 500 750 1000 1250)

home=/home1/zhuyixin/sc2/ImmAssm

while getopts 's:p:a:m:h:' option
do
    case "$option" in
        s) species=${OPTARG};;
        p) chr1=${OPTARG};;
        a) chr2=${OPTARG};;
        m) chr3=${OPTARG};;
        h) chr4=${OPTARG};;
    esac
done

# Run the processing script with the appropriate threshold
bash ${home}/code/cigar.sh ${THRESHOLDS[$SLURM_ARRAY_TASK_ID]} $species $chr1 $chr2 $chr3 $chr4