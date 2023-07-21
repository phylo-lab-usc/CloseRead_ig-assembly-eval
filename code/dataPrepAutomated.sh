#!/bin/bash
#SBATCH --job-name=automated    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log


module load python
module load conda
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly
source /etc/profile.d/modules.sh
module load gcc/11.3.0
module load samtools/1.17

cat code/name.txt | while read line
do 
    #check if this species' data is in bam or not, convert to fastq if yes
    count=`ls -1 hifi_fastq/${line}/*.bam 2>/dev/null | wc -l`
    if [ $count != 0 ]
    then 
    for f in hifi_fastq/${line}/*.bam ; do 
	(
	name=${f%.ccs*};
	echo $name;
       	bam2fastq -o $name -j 5 $f ;
	)&
    done ; wait
    fi
    #merge all the fastq file into a single file
    cat hifi_fastq/${line}/*.fastq > hifi_fastq/${line}/${line}_merged.fastq
    #map the merged fastq file to the coresponding assembly
    minimap2 -t 20 -a /home1/zhuyixin/sc1/ImmAssm/assemblies/${line}*.fasta /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/${line}/${line}_merged.fastq > /home1/zhuyixin/sc1/ImmAssm/aligned_sam/${line}/${line}_merged.sam
    #convert the SAM result to sorted BAM format
    samtools sort /home1/zhuyixin/sc1/ImmAssm/aligned_sam/${line}/${line}_merged.sam -o /home1/zhuyixin/sc1/ImmAssm/aligned_bam/${line}/${line}_merged_sorted.bam
    #index the sorted BAM file
    samtools index /home1/zhuyixin/sc1/ImmAssm/aligned_bam/${line}/${line}_merged_sorted.bam
done

 



