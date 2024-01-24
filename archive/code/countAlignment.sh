#!/bin/bash
#SBATCH --job-name=countAlignment    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/countAlignment%j.log   # Standard output and error log
#SBATCH --mem=200G



source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly


while getopts 'g:e:' option
do
    case "$option" in        
        g) genome=${OPTARG};;
        e) extended=${OPTARG};;
    esac
done

HOME=/home1/zhuyixin/sc1/ImmAssm
cat ${HOME}/code/name.txt | while read line
do 
    echo ${line}
    #awk -v OFS='\t' {'print $1,$2'} assemblies/${line}.merged.fasta.fai > assemblies/${line}.genome
    echo "finished genome prep"
    #extract secondary alignment
    #samtools view -b -f 256 -@ 40 aligned_bam/${genome}/${line}/${line}_merged*.bam > split_bam/${genome}/${line}_merged_secondary.bam
    echo "finished secondary alignment"
    #extract primary alignment 
    #samtools view -b -F 0x800 -F 0x100 -@ 40 aligned_bam/${genome}/${line}/${line}_merged*.bam > split_bam/${genome}/${line}_merged_primary.bam
    echo "finished primary alignment"
    #bedtools to generate stats
    #coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/functional/${line}/gene_IGH_pos_sorted.bed -b split_bam/${genome}/${line}_merged_primary.bam > alignment_count/${extended}/${genome}/${line}_funcional_primary_count.txt
    coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/functional/${line}_functional.bed -b split_bam/${genome}/${line}_merged_primary.bam > alignment_count/${extended}/${genome}/${line}_funcional_primary_count.txt
    echo "1"
    #coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/functional/${line}/gene_IGH_pos_sorted.bed -b split_bam/${genome}/${line}_merged_secondary.bam > alignment_count/${extended}/${genome}/${line}_funcional_secondary_count.txt 
    coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/functional/${line}_functional.bed -b split_bam/${genome}/${line}_merged_secondary.bam > alignment_count/${extended}/${genome}/${line}_funcional_secondary_count.txt
    echo "2"
    #coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/nonfunctional/${line}/gene_IGH_pos_sorted.bed -b split_bam/${genome}/${line}_merged_primary.bam > alignment_count/${extended}/${genome}/${line}_nonfuncional_primary_count.txt
    coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/nonfunctional/${line}_nonfunctional.bed -b split_bam/${genome}/${line}_merged_primary.bam > alignment_count/${extended}/${genome}/${line}_nonfuncional_primary_count.txt
    echo "3"
    #coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/nonfunctional/${line}/gene_IGH_pos_sorted.bed -b split_bam/${genome}/${line}_merged_secondary.bam > alignment_count/${extended}/${genome}/${line}_nonfuncional_secondary_count.txt
    coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/nonfunctional/${line}_nonfunctional.bed -b split_bam/${genome}/${line}_merged_secondary.bam > alignment_count/${extended}/${genome}/${line}_nonfuncional_secondary_count.txt  
    echo "4"
done

