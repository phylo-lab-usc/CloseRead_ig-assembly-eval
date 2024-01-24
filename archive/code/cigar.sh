#!/bin/bash
#SBATCH --job-name=cigar    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60                    # Run on a single CPU 
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/hcigar_%j.log   # Standard output and error log
#SBATCH --mem=100G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

home=/home1/zhuyixin/sc2/ImmAssm
while getopts 's:c:' option
do
    case "$option" in
        s) species=${OPTARG};;
        c1) chr1=${OPTARG};;
        c2) chr2=${OPTARG};;
    esac
done

samtools view -@ 60 -b ${home}/aligned_bam/${species}_merged_sorted.bam ${c1} ${c2} > ${home}/aligned_bam/${species}.ig.bam
samtools view -F 0x800 -F 0x100 -@ 60 ${home}/aligned_bam/${species}.ig.bam > ${home}/aligned_bam/${species}.ig.pri.sam

for thre in 500 750 1000 1250; do
    cat ${home}/aligned_bam/${species}.ig.pri.sam | while read readID readFlag chr pos mapq cigar rnext pnext tlen seq qual tagMS tag
    do
        long_blocks_sum=0
        long_blocks_sum=$(echo $cigar | grep -oE '[0-9]+M' | awk '{ split($0, a, "M"); if(a[1] >= ${thre}) sum += a[1] } END { print sum+0 }')
        total_length=$(echo $cigar | grep -oE '[0-9]+[MDI]' | awk '{ split($0, a, /[MDI]/); sum += a[1] } END { print sum }')
        fraction=0
        if [ $total_length -ne 0 ]; then
            fraction=$(awk "BEGIN {print $long_blocks_sum/$total_length}")
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${readID}" "${chr}" "${pos}" "${mapq}" "${total_length}" "${fraction}" >> ${home}/cigar_result/${species}.${chr}.pri.long_blocks.${thre}.txt
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${readID}" "${chr}" "${pos}" "${mapq}" "${total_length}" "${fraction}" >> ${home}/cigar_result/${species}.ig.pri.long_blocks.${thre}.txt
    done
    python3 ${home}/code/analyze_distributions.py "${home}/cigar_result/${species}.${c1}.pri.long_blocks.${thre}.txt" "$thre" "${species}" "${c1}"
    python3 ${home}/code/analyze_distributions.py "${home}/cigar_result/${species}.${c2}.pri.long_blocks.${thre}.txt" "$thre" "${species}" "${c2}"
    python3 ${home}/code/analyze_distributions.py "${home}/cigar_result/${species}.ig.pri.long_blocks.${thre}.txt" "$thre" "${species}" "ig"

done



