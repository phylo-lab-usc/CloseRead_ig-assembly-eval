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

species=mCynVol1
HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality
outdir=mCynVol1.corrected
ref=${HOME}/assemblies/${species}.corr.merged.fasta
fastq=${HOME}/hifi_fastq/${species}/${species}_merged.fastq

minimap2 -t 32 -ax map-hifi $ref $fastq > ${HOME}/aligned_sam/${outdir}/${species}_merged.sam
#convert the SAM result to sorted BAM format
echo "converting SAM to sorted BAM"
samtools sort -@ 32 ${HOME}/aligned_sam/${outdir}/${species}_merged.sam -o ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam
#index the sorted BAM file
echo "indexing sorted BAM"
samtools index -c -@ 32 ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam
samtools view -b -F 0x800 -F 0x100 -@ 32 ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam > ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted_primary.bam
samtools index -c -@ 32 ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted_primary.bam
