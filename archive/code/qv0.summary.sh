#!/bin/bash
#SBATCH --job-name=qv0.summary    # Job name
#SBATCH --nodes=1
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/qv0.summary_%j.log   # Standard output and error log
#SBATCH --mem=50G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

cd /home1/zhuyixin/sc2/ImmAssm/aligned_bam

#samtools view -F 0x800 -F 0x100 -@ 60 rheMac_merged_sorted.bam > rheMac_merged_sorted.primary.sam
#samtools view -f 256 -@ 60 rheMac_merged_sorted.bam > rheMac_merged_sorted.secondary.sam

#check if # of reads match # of lines greped from sam file - can do additional unique check
#LC_ALL=C grep -F -f ${species}.unique.qv0.readID.txt ${species}_merged_sorted.primary.sam > ${species}.qv0.pri.sam
#LC_ALL=C grep -F -f ${species}.unique.qv0.readID.txt ${species}_merged_sorted.secondary.sam > ${species}.qv0.sec.sam

cat ${species}.unique.qv0.readID.txt | while read readID
do
    date
    fname=$(echo $readID | tr / :)
    LC_ALL=C grep -F $readID ${species}.qv0.pri.sam > ${species}/${fname}.primary.sam
    date
    temp=$(awk '{ print $21 }' ${species}/${fname}.primary.sam)
    nSeed=${temp#rl:i:}
    chr=$(awk '{ print $3 }' ${species}/${fname}.primary.sam)
    start=$(awk '{ print $4 }' ${species}/${fname}.primary.sam)
    echo $nSeed
    echo $chr
    echo $start
    LC_ALL=C grep -F $readID ${species}.qv0.sec.sam > ${species}/${fname}.secondary.sam
    date
    nSec=$(wc -l < "${species}/${fname}.secondary.sam")
    printf "%s\t%s\t%s\t%s\t%s\t" "${readID}" "${chr}" "${start}" "${nSeed}" "${nSec}" >> ${species}.qv0.summary.txt
    chrSec=($(awk '{ print $3 }' ${species}/${fname}.secondary.sam))
    startSec=($(awk '{ print $4 }' ${species}/${fname}.secondary.sam))
    for i in $( eval echo {0..$((nSec - 1))})
    do
        printf "%s\t%s\t" "${chrSec[i]}" "${startSec[i]}" >> ${species}.qv0.summary.txt
    done
    printf "\n" >> ${species}.qv0.summary.txt
done



