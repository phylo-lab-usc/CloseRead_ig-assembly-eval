#!/bin/bash
#SBATCH --job-name=minimap2    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/minimap2_%j.log   # Standard output and error log
#SBATCH --mem=200G

source /etc/profile.d/modules.sh
module load gcc/11.3.0
module load samtools/1.17
module load python
module load conda
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly
cd ~

home=/home1/zhuyixin/sc1/ImmAssm
cat ${home}/gene_position/IgH_loci.txt | while read -r name chrom start end
do
    samtools faidx ${home}/assemblies/${name}.merged.fasta ${chrom}:${start}-${end} > ${home}/assemblies/${chrom}_${name}.txt
done