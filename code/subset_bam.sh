#!/bin/bash
#SBATCH --job-name=subset_bam    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/subset_bam_%j.log   # Standard output and error log
#SBATCH --mem=200G



source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

while getopts 'e:s:f:' option
do
    case "$option" in
        e) extended=${OPTARG};;
        s) species=${OPTARG};;
        f) function=${OPTARG};;
    esac
done

HOME=/home1/zhuyixin/sc2/ImmAssm
BED=${HOME}/gene_position/${extended}/${species}/${species}_${function}.bed

cat $BED | while read chrom start end sv score strand
do
    samtools view -b -@ 60 ${HOME}/aligned_bam/${species}_merged_sorted.bam "${chrom}:${start}-${end}" > ${HOME}/aligned_bam/${species}/${species}_${chrom}:${start}-${end}.bam
    samtools view -F 0x800 -F 0x100 -@ 60 ${HOME}/aligned_bam/${species}/${species}_${chrom}:${start}-${end}.bam | awk '$5==0 { print $0 }' > ${HOME}/aligned_bam/${species}/${species}_${chrom}:${start}-${end}.qv0.bam
    nRow=$(wc -l < "${HOME}/aligned_bam/${species}/${species}_${chrom}:${start}-${end}.qv0.bam")
    if [ $nRow != 0 ]
    then
        cat ${HOME}/aligned_bam/${species}/${species}_${chrom}:${start}-${end}.qv0.bam | awk '{ print $1 }' >> ${HOME}/aligned_bam/${species}.qv0.readID.txt
    fi
done

cat ${HOME}/aligned_bam/${species}.qv0.readID.txt | sort | uniq > ${HOME}/aligned_bam/${species}.unique.qv0.readID.txt 

