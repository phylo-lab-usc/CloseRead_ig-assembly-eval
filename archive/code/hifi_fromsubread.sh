#!/bin/bash
#SBATCH --job-name=hifi_fromsubread    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/hifi_fromsubread_%j.log   # Standard output and error log
#SBATCH --mem=100G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

home=/home1/zhuyixin/sc2/ImmAssm

cat ${home}/code/name.txt | while read line
do
    for f in ${home}/sim_sam/${line}/*6.sam
    do
        name=${f%.sam};
        samtools view -bS -@ 40 ${f} > ${name}.bam
        ccs -j 40 ${name}.bam ${name}.fastq.gz
    done
done