#!/bin/bash
#SBATCH --job-name=qv0.summary    # Job name
#SBATCH --nodes=1
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/qv0.summary_%j.log   # Standard output and error log
#SBATCH --mem=20G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

cd /home1/zhuyixin/sc2/ImmAssm/aligned_bam

while getopts 'n:' option
do
    case "$option" in
        n) iter=${OPTARG};;
    esac
done

cat human_ID/human.14.ID.${iter} | while read readID
do
    date
    fname=$(echo $readID | tr / :)
    LC_ALL=C grep -F $readID human.chr14.qv0.pri.sam > human/${fname}.chr14.primary.sam
    date
    temp=$(awk '{ print $21 }' human/${fname}.chr14.primary.sam)
    nSeed=${temp#rl:i:}
    chr=$(awk '{ print $3 }' human/${fname}.chr14.primary.sam)
    start=$(awk '{ print $4 }' human/${fname}.chr14.primary.sam)
    echo $nSeed
    echo $chr
    echo $start
    LC_ALL=C grep -F $readID human.chr14.qv0.sec.sam > human/${fname}.chr14.secondary.sam
    date
    nSec=$(wc -l < "human/${fname}.chr14.secondary.sam")
    printf "%s\t%s\t%s\t%s\t%s\t" "${readID}" "${chr}" "${start}" "${nSeed}" "${nSec}" >> human.14.qv0.summary.${iter}.txt
    chrSec=($(awk '{ print $3 }' human/${fname}.chr14.secondary.sam))
    startSec=($(awk '{ print $4 }' human/${fname}.chr14.secondary.sam))
    for i in $( eval echo {0..$((nSec - 1))})
    do 
        printf "%s\t%s\t" "${chrSec[i]}" "${startSec[i]}" >> human.14.qv0.summary.${iter}.txt
    done
    printf "\n" >> human.14.qv0.summary.${iter}.txt
done 