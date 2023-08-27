#!/bin/bash
#SBATCH --job-name=merge_automated    # Job name
#SBATCH --nodes=1
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=merge_automated.log   # Standard output and error log


#HOME=/home1/zhuyixin/sc1/ImmAssm
while getopts 'f:p:e:' option
do
    case "$option" in
        f) function=${OPTARG};;
        p) HOME=${OPTARG};;
        e) extended=${OPTARG};;
    esac
done
mkdir ${HOME}/gene_position/${extended}/combined
mkdir ${HOME}/gene_position/${extended}/combined/${function}
cat ${HOME}/code/name.txt | while read line
do 
    cat ${HOME}/gene_position/${extended}/primary/${function}/${line}/gene_IGH_pos_sorted.bed ${HOME}/gene_position/${extended}/alt/${function}/${line}/gene_IGH_pos_sorted.bed > ${HOME}/gene_position/${extended}/combined/${function}/${line}_${function}.bed
done