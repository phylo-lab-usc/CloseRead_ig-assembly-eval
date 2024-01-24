#!/bin/bash
#SBATCH --job-name=minimap2    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/minimap2_%j.log   # Standard output and error log
#SBATCH --mem=200G

module load python
module load conda
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly
cd ~

home=/home1/zhuyixin/sc1/ImmAssm
for f in ~/sc1/ImmAssm/assemblies/*.alt.fasta ;
do
    temp=${f#*assemblies/} 
    name=${temp%.alt*};
    seqtk seq -F '#' ${home}/assemblies/${name}.alt.fasta > ${home}/assemblies/${name}.alt.fastq
    minimap2 -t 40 ${home}/assemblies/${name}.fasta ${home}/assemblies/${name}.alt.fastq > ${home}/chrom_match/${name}.paf
done