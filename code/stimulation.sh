#!/bin/bash
#SBATCH --job-name=stimulation    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40                    # Run on a single CPU 
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=log/stimulation_%j.log   # Standard output and error log
#SBATCH --mem=200G

source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

home=/home1/zhuyixin/sc2/ImmAssm
cat ${home}/code/name.txt | while read line
do
    mkdir ${home}/sim_sam/${line}
    cd ${home}/sim_sam/${line}
    ${home}/code/pbsim --strategy wgs --method qshmm --qshmm ~/bin/pbsim3/data/QSHMM-RSII.model --difference-ratio 22:45:33 --depth 30 --seed 10010101 --genome ${home}/assemblies/${line}.merged.fasta --pass-num 15
done
