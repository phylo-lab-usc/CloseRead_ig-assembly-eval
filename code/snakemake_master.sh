#!/bin/bash
#SBATCH --job-name=snakemake_master    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10                    # Run on a single CPU 
#SBATCH --time=12:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/snakemake_master%j.log   # Standard output and error log
#SBATCH --mem=100G


source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

# Run Snakemake (submitting jobs to SLURM)
snakemake --cluster "sbatch -A mpennell_978 -p gpu --ntasks=1 --cpus-per-task=32 --output=log/%j.out --time=24:00:00 --mem=65GB" \
          --snakefile Snakefile \
          --printshellcmds --reason --verbose --latency-wait 60000 --cores all --jobs 2 --rerun-incomplete
