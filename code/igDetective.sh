#!/bin/bash
#SBATCH --job-name=igDetect    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=igDetective%j.log   # Standard output and error log
#SBATCH --mem=60G




source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/IGdetective
conda env list

cd /home1/zhuyixin/IgDetective
/home1/zhuyixin/.conda/envs/IGdetective/bin/python run_iterative_igdetective.py ~/sc1/ImmAssm/assemblies/mMelMel3.fasta ~/sc1/ImmAssm/mammalian_igdetective_v2.0/mMelMel3_igdetective/