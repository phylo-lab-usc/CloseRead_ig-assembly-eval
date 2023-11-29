#!/bin/bash
#SBATCH --job-name=stimulation    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60                    # Run on a single CPU 
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/stimulation_%j.log   # Standard output and error log
#SBATCH --mem=100G

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
#mkdir ${home}/test/${species}
#cd ${home}/test/${species}
cd ${home}/sim_sam/${species}
${home}/code/pbsim --strategy wgs --method qshmm --qshmm ~/bin/pbsim3/data/QSHMM-RSII.model --difference-ratio 22:45:33 --depth 30 --genome ${home}/assemblies/${species}/NC_000014.9.fasta --prefix 14 --pass-num 10
now=$(date)
echo "$now"

<<COMMENT1
for f in ${home}/test/${line}/*.sam
do
    name=${f%.sam};
    samtools view -bS -@ 40 ${f} > ${name}.bam
    rm -rf ${f}
    ccs -j 40 ${name}.bam ${name}.fastq.gz
    rm -rf ${name}.bam
done
COMMENT1

samtools view -bS -@ 60 14*sam > 14.bam
rm -rf 14*sam
ccs -j 60 14.bam 14.fastq.gz

now=$(date)
echo "$now"
