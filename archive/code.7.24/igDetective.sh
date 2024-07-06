#!/bin/bash
#SBATCH --job-name=igDetect    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/igDetective%j.log   # Standard output and error log
#SBATCH --mem=30G

source /etc/profile.d/modules.sh
module load gcc/11.3.0
module load samtools/1.17
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda config --append envs_dirs /home1/zhuyixin/.conda/envs/
conda activate /home1/zhuyixin/.conda/envs/IGdetective
conda env list
conda info --envs
echo "PATH: $PATH"
which python

/home1/zhuyixin/.conda/envs/IGdetective/bin/python3 /home1/zhuyixin/IgDetective/run_iterative_igdetective.py $1 $2
touch /home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/${3}.${4}.txt
