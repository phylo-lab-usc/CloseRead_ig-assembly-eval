#!/bin/bash
#SBATCH --job-name=igv    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=igv.log   # Standard output and error log

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

function igv {
    mkdir ${HOME}/igv/${genome}_${extended}
    mkdir ${HOME}/igv/${genome}_${extended}/${function}
    mkdir ${HOME}/igv/${genome}_${extended}/${function}/${1}
    cat ${HOME}/gene_position/${extended}/${genome}/${function}/${1}_${function}.bed | while read chrom start end sv score strand
    do
    echo $sv
    mkdir ${HOME}/igv/${genome}_${extended}/${function}/${1}/snapshot
    bam="${HOME}/aligned_bam/${genome}/${1}/${1}_merged_sorted.bam"
    echo "new"
    echo "preference SAM.SHOW_ALL_BASES 0"
    echo "preference SAM.SHOW_SOFT_CLIPPED true"
    echo "sort QUALITY"
    if [ $genome == 'combined' ]
    then
        echo "genome ${HOME}/assemblies/${1}.merged.fasta"
    else   
        echo "genome ${HOME}/assemblies/${1}.fasta"
    fi
    echo "snapshotDirectory ${HOME}/igv/${genome}_${extended}/${function}/${1}/snapshot"
    echo "load ${bam}"
    echo "load ${HOME}/gene_position/${extended}/${genome}/${function}/${1}_${function}.bed"
    echo "goto ${chrom}:${start}-${end}"
    echo "snapshot ${1}_merged_${sv}.png"
    done    > /home1/zhuyixin/sc1/ImmAssm/igv/${genome}_${extended}/${function}/${1}/igv.txt
}

HOME=/home1/zhuyixin/sc1/ImmAssm
while getopts 'f:g:e:' option
do
    case "$option" in
        f) function=${OPTARG};;
        g) genome=${OPTARG};;
        e) extended=${OPTARG};;
    esac
done
cat ${HOME}/code/name.txt | while read line
do 
    igv $line $function $genome $extended
done
