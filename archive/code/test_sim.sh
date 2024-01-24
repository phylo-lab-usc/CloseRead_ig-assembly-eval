#!/bin/bash
#SBATCH --job-name=test_sim    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60                    # Run on a single CPU 
#SBATCH --time=120:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/test_sim_%j.log   # Standard output and error log
#SBATCH --mem=60G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

now=$(date)
echo "$now"
home=/home1/zhuyixin/sc2/ImmAssm
while getopts s: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
    esac
done
echo ${species}
mkdir ${home}/test/${species}
cd ${home}/test/${species}
for i in {7..10}
do 
    for a in ${home}/assemblies/${species}/*fasta
    do 
        temp=${a%.fasta}
        name=${temp#*${species}/}
        echo $name
        ${home}/code/pbsim --strategy wgs --method qshmm --qshmm ~/bin/pbsim3/data/QSHMM-RSII.model --difference-ratio 22:45:33 --depth 3 --genome ${a} --pass-num 10 --prefix ${i}_${name} --seed $i --id-prefix ${i}_${name}
        for f in ${home}/test/${species}/${i}_*.sam
        do
            name1=${f%.sam}
            samtools view -bS -@ 60 ${f} > ${name1}.bam
            rm -rf $f
            rm -rf ${name1}.maf
            rm -rf ${name1}.ref
            ccs -j 60 ${name1}.bam ${name1}.fastq.gz
            rm -rf ${name1}.bam
        done
    done
done


now=$(date)
echo "$now"
