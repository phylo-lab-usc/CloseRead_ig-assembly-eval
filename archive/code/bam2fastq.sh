#!/bin/bash
#SBATCH --job-name=bam2fastq    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

module load python
module load conda
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly
cd ~/sc1/ImmAssm/hifi_bam/mCanLor1
N=5
for f in *.bam ; do 
	(
	 name=${f%.ccs*};
	 echo $name;
       	 bam2fastq -o $name -j 5 $f ;
	)&
done ; wait
