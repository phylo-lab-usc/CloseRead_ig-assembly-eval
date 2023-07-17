#!/bin/bash
#SBATCH --job-name=minimap2    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

module load python
module load conda
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly
cd ~
minimap2 -t 20 -a /home1/zhuyixin/sc1/ImmAssm/assemblies/mCanLor1.pri.cur.20210315.fasta /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mCanLor1/mCanLor1_merged.fastq > /home1/zhuyixin/sc1/ImmAssm/aligned_sam/mCanLor1/mCanLor1_merged.sam