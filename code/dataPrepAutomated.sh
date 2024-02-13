#!/bin/bash
#SBATCH --job-name=automated    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/dataPrep%j.log   # Standard output and error log
#SBATCH --mem=100G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/23.10.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

while getopts s:c:w: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        c) correction=${OPTARG};;
        w) source=${OPTARG};;
    esac
done
HOME=/home1/zhuyixin/sc1/AssmQuality

#check if this species' data is in bam or not, convert to fastq if yes
count=`ls -1 ${HOME}/${source}/${species}/*.bam 2>/dev/null | wc -l`
merged=`ls -1 ${HOME}/${source}/${species}/*_merged.fastq 2>/dev/null | wc -l`
if [ $count != 0 ] && [ $merged == 0 ]
then 
    echo "bam file exists, converting to fastq"
    for f in ${HOME}/${source}/${species}/*.bam ; do 
    (
        name=${f%.ccs*};
        echo $f;
        echo $name;
        bam2fastq -o $name -j 60 $f ;
    )&
    done ; wait
    echo "bam to fastq conversion done"
fi
if [ $merged == 0 ]
then
    #unzip fastq gz file
    echo "gunziping fastq files"
    gunzip ${HOME}/${source}/${species}/*.gz
    #merge all the fastq file into a single file
    echo "merging fastq files"
    cat ${HOME}/${source}/${species}/*.fastq > ${HOME}/${source}/${species}/${species}_merged.fastq
    cat ${HOME}/${source}/${species}/*.fq >> ${HOME}/${source}/${species}/${species}_merged.fastq
fi
#create output directories
mkdir ${HOME}/aligned_sam/
mkdir ${HOME}/aligned_bam/
if [ $correction == "True" ]
then
    mkdir ${HOME}/aligned_sam/${species}.corrected
    mkdir ${HOME}/aligned_bam/${species}.corrected
    #map the merged fastq file to the coresponding assembly
    echo "mapping fastq to assembly"
    minimap2 -t 60 -a ${HOME}/assemblies/${species}.ljacorr.merged.fasta ${HOME}/${source}/${species}/${species}_merged.fastq > ${HOME}/aligned_sam/${species}.corrected/${species}_merged.sam
    #convert the SAM result to sorted BAM format
    echo "converting SAM to sorted BAM"
    samtools sort -@ 60 ${HOME}/aligned_sam/${species}.corrected/${species}_merged.sam -o ${HOME}/aligned_bam/${species}.corrected/${species}_merged_sorted.bam
    #index the sorted BAM file
    echo "indexing sorted BAM"
    samtools index ${HOME}/aligned_bam/${species}.corrected/${species}_merged_sorted.bam
else
    mkdir ${HOME}/aligned_sam/${species}
    mkdir ${HOME}/aligned_bam/${species}
    #map the merged fastq file to the coresponding assembly
    echo "mapping fastq to assembly"
    minimap2 -t 60 -a ${HOME}/assemblies/${species}.merged.fasta ${HOME}/${source}/${species}/${species}_merged.fastq > ${HOME}/aligned_sam/${species}/${species}_merged.sam
    #convert the SAM result to sorted BAM format
    echo "converting SAM to sorted BAM"
    samtools sort -@ 60 ${HOME}/aligned_sam/${species}/${species}_merged.sam -o ${HOME}/aligned_bam/${species}/${species}_merged_sorted.bam
    #index the sorted BAM file
    echo "indexing sorted BAM"
    samtools index ${HOME}/aligned_bam/${species}/${species}_merged_sorted.bam
fi






