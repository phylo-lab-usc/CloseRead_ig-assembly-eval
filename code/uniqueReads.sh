#!/bin/bash
#SBATCH --job-name=unique60    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/unique_read%j.log   # Standard output and error log
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

for f in ~/sc1/ImmAssm/split_bam/${genome}/*_merged_primary.bam;
do
    temp=${f#*${genome}/}
    name=${temp%_merged*}
    echo $name
    if [ "$genome" = "combined" ]; then
        samtools view -b -F 256 -q 60 -@ 40 $f > ~/sc1/ImmAssm/split_bam/combined/${name}_unique60_primary.bam
        coverageBed -counts -sorted -nobuf -g assemblies/${name}.genome -a gene_position/${extended}/${genome}/functional/${name}_functional.bed -b split_bam/${genome}/${name}_unique60_primary.bam > alignment_count/${extended}/${genome}/${line}_funcional_unique60_primary_count.txt
        coverageBed -counts -sorted -nobuf -g assemblies/${name}.genome -a gene_position/${extended}/${genome}/nonfunctional/${name}_nonfunctional.bed -b split_bam/${genome}/${name}_unique60_primary.bam > alignment_count/${extended}/${genome}/${line}_nonfuncional_unique60_primary_count.txt
    else
        coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/functional/${line}/gene_IGH_pos_sorted.bed -b split_bam/${genome}/${name}_unique60_primary.bam > alignment_count/${extended}/${genome}/${line}_funcional_unique60_primary_count.txt
        coverageBed -counts -sorted -nobuf -g assemblies/${line}.genome -a gene_position/${extended}/${genome}/nonfunctional/${line}/gene_IGH_pos_sorted.bed -b split_bam/${genome}/${name}_unique60_primary.bam > alignment_count/${extended}/${genome}/${line}_nonfuncional_unique60_primary_count.txt
    fi
done
