#!/bin/bash
#SBATCH --job-name=automated    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/dataPrep%j.log   # Standard output and error log
#SBATCH --mem=60G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

while getopts s:w:h: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        w) source=${OPTARG};;
        h) haploid=${OPTARG};;
    esac
done
HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality

echo $haploid