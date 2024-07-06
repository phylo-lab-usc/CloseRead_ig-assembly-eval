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

while getopts s:w:h:d: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        w) source=${OPTARG};;
        h) haploid=${OPTARG};;
        d) HOME=${OPTARG};;
    esac
done

#check if this species' raw data has .bam format or not, convert to fastq if yes
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
        bam2fastq -o $name -j 32 $f ;
    )&
    done ; wait
    echo "bam to fastq conversion done"
fi


#If .merged.fastq exists, skip. Otherwise, gunzip fastq files and concatenate to .merged.fastq
if [ $merged == 0 ]
then
    #unzip fastq gz file
    echo "gunziping fastq files"
    gunzip -f ${HOME}/${source}/${species}/*.gz
    #merge all the fastq file into a single file
    echo "merging fastq files"
    cat ${HOME}/${source}/${species}/*.fastq > ${HOME}/${source}/${species}/${species}_merged.fastq
    cat ${HOME}/${source}/${species}/*.fq >> ${HOME}/${source}/${species}/${species}_merged.fastq
fi


#define output and input
outdir=${species}
if [ "$haploid" == "True" ]
then
    ref=${HOME}/assemblies/${species}.pri.fasta
else
    #ref=${HOME}/assemblies/${species}.merged.fasta
    ref=${HOME}/assemblies/${species}.merged.fasta
fi

#create output directories
mkdir ${HOME}/aligned_sam/
mkdir ${HOME}/aligned_bam/
mkdir ${HOME}/aligned_sam/${outdir}
mkdir ${HOME}/aligned_bam/${outdir}
fastq=${HOME}/${source}/${species}/${species}_merged.fastq

#map the merged fastq file to the coresponding assembly
echo "mapping fastq to assembly"
minimap2 -t 32 -ax map-hifi $ref $fastq > ${HOME}/aligned_sam/${outdir}/${species}_merged.sam
#convert the SAM result to sorted BAM format
echo "converting SAM to sorted BAM"
samtools sort -@ 32 ${HOME}/aligned_sam/${outdir}/${species}_merged.sam -o ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam
#index the sorted BAM file
echo "indexing sorted BAM"
samtools index -c -@ 32 ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam
#rm -rf ${HOME}/aligned_sam/${outdir}/${species}_merged.sam




