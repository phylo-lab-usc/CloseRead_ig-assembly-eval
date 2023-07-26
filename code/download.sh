#!/bin/bash
#SBATCH --job-name=download    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=download.log   # Standard output and error log

#module load python
#module spider anaconda
#module load conda
#conda init
#source /spack/apps/anaconda3/2021.05/etc/profile.d/conda.sh
#conda activate igortest

mkdir /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mMacEug1
cd /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mMacEug1

wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m54306U_210508_040324.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m54306U_210509_130224.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m54306Ue_210702_174342.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m54306Ue_210704_043727.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m54306Ue_211001_163112.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m64330e_211013_010414.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m64334e_211013_013414.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Macropus_eugenii/mMacEug1/genomic_data/pacbio_hifi/m64334e_211027_184708.hifi_reads.fastq.gz
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201003_145620.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201004_152428.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201014_175437.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201015_180930.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi

