#!/bin/bash
#SBATCH --job-name=automated    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=dataPrep.log   # Standard output and error log




source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

while getopts g: flag
do
    case "${flag}" in
        g) genome=${OPTARG};;
    esac
done
echo $genome
HOME=/home1/zhuyixin/sc1/ImmAssm
cat ${HOME}/code/name.txt | while read line
do 
    #check if this species' data is in bam or not, convert to fastq if yes
    count=`ls -1 ${HOME}/hifi_fastq/${line}/*.bam 2>/dev/null | wc -l`
    merged=`ls -1 ${HOME}/hifi_fastq/${line}/*_merged.fastq 2>/dev/null | wc -l`
    if [ $count != 0 ] && [ $merged == 0 ]
    then 
        for f in ${HOME}/hifi_fastq/${line}/*.bam ; do 
        (
            name=${f%.ccs*};
            echo $f;
            echo $name;
            bam2fastq -o $name -j 5 $f ;
        )&
        done ; wait
    fi
    if [ $merged == 0 ]
    then
        #unzip fastq gz file
        gunzip ${HOME}/hifi_fastq/${line}/*.gz
        #merge all the fastq file into a single file
        cat ${HOME}/hifi_fastq/${line}/*.fastq > ${HOME}/hifi_fastq/${line}/${line}_merged.fastq
    fi
    #create output directories
    mkdir ${HOME}/aligned_sam/${genome}/${line}
    mkdir ${HOME}/aligned_bam/${genome}/${line}
    #map the merged fastq file to the coresponding assembly
    minimap2 -t 20 -a ${HOME}/assemblies/${line}*merged.fasta ${HOME}/hifi_fastq/${line}/${line}_merged.fastq > ${HOME}/aligned_sam/${genome}/${line}/${line}_merged.sam
    #convert the SAM result to sorted BAM format
    samtools sort ${HOME}/aligned_sam/${genome}/${line}/${line}_merged.sam -o ${HOME}/aligned_bam/${genome}/${line}/${line}_merged_sorted.bam
    #index the sorted BAM file
    samtools index ${HOME}/aligned_bam/${genome}/${line}/${line}_merged_sorted.bam
done

 



