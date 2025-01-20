#!/bin/bash
while getopts s:w:h:d:t: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        w) source=${OPTARG};;
        h) haploid=${OPTARG};;
        d) HOME=${OPTARG};;
        t) threads=${OPTARG};;
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
        name=${f%.bam*};
        echo $f;
        echo $name;
        pbindex $f;
        bam2fastq -o $name -j $threads $f ;
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
minimap2 -t $threads -ax map-hifi $ref $fastq > ${HOME}/aligned_sam/${outdir}/${species}_merged.sam
#convert the SAM result to sorted BAM format
echo "converting SAM to sorted BAM"
samtools sort -@ $threads ${HOME}/aligned_sam/${outdir}/${species}_merged.sam -o ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam
#index the sorted BAM file
echo "indexing sorted BAM"
samtools index -c -@ $threads ${HOME}/aligned_bam/${outdir}/${species}_merged_sorted.bam
#rm -rf ${HOME}/aligned_sam/${outdir}/${species}_merged.sam




