#!/bin/bash
#SBATCH --job-name=coverage    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30                    # Run on a single CPU 
#SBATCH --time=6:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/coverage%j.log   # Standard output and error log
#SBATCH --mem=100G


source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly



while getopts s:a:b: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        a) assemblies=${OPTARG};;
        b) bam=${OPTARG};;
    esac
done
HOME=/home1/zhuyixin/sc1/AssmQuality

echo ${bam}
for gene in IGH IGK IGL ; do
    echo $gene
    most_freq_chr_pri=$(awk '{print $1}' ${HOME}/gene_position/${species}/pri/${species}_${gene}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 1 | awk '{print $2}')
    most_freq_chr_alt=$(awk '{print $1}' ${HOME}/gene_position/${species}/alt/${species}_${gene}_pos.sorted.bed | sort | uniq -c | sort -nr | head -n 1 | awk '{print $2}')
    priloc=$(awk -v chr="$most_freq_chr_pri" 'BEGIN{min=""; max=0} $1 == chr {if (min=="" || $2<min) min=$2; if ($3>max) max=$3} END {print chr":"min"-"max}' ${HOME}/gene_position/${species}/pri/${species}_${gene}_pos.sorted.bed)
    altloc=$(awk -v chr="$most_freq_chr_alt" 'BEGIN{min=""; max=0} $1 == chr {if (min=="" || $2<min) min=$2; if ($3>max) max=$3} END {print chr":"min"-"max}' ${HOME}/gene_position/${species}/alt/${species}_${gene}_pos.sorted.bed)
    echo $priloc $altloc
    samtools mpileup -Q 0 -q 0 -aa -f ${assemblies} -r $priloc ${bam} > ${HOME}/errorStats/${species}/${gene}_pri_pileup.txt
    samtools mpileup -Q 0 -q 0 -aa -f ${assemblies} -r $altloc ${bam} > ${HOME}/errorStats/${species}/${gene}_alt_pileup.txt
    #/home1/zhuyixin/.conda/envs/assembly/bin/python ${HOME}/code/coverageAnalysis.py ${HOME}/errorStats/${species}/${gene}_pri_pileup.txt ${HOME}/errorStats/${species}/${gene}_pri_coverage.txt ${species} ${gene}_pri
    #/home1/zhuyixin/.conda/envs/assembly/bin/python ${HOME}/code/coverageAnalysis.py ${HOME}/errorStats/${species}/${gene}_alt_pileup.txt ${HOME}/errorStats/${species}/${gene}_alt_coverage.txt ${species} ${gene}_alt
done