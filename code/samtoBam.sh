#!/bin/bash
#SBATCH --job-name=samtobam    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=samtobam_%j.log   # Standard output and error log

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

while getopts s: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
    esac
done
samtools sort -@ 20 ~/sc1/ImmAssm/aligned_sam/primary/${species}/${species}*_merged.sam -o ~/sc1/ImmAssm/aligned_bam/primary/${species}/${species}_merged.bam
samtools index ~/sc1/ImmAssm/aligned_bam/primary/${species}/${species}_merged.bam
