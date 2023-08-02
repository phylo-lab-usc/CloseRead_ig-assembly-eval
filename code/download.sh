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

cd /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mNeoNeb1
wget https://genomeark.s3.amazonaws.com/species/Neofelis_nebulosa/mNeoNeb1/genomic_data/pacbio_hifi/m64055e_211110_175241.hifi_reads.fastq.gz

mkdir /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mCynVol1
cd /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mCynVol1

wget https://genomeark.s3.amazonaws.com/species/Cynocephalus_volans/mCynVol1/genomic_data/pacbio_hifi/m54306Ue_211202_201212.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Cynocephalus_volans/mCynVol1/genomic_data/pacbio_hifi/m54306Ue_211213_200140.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Cynocephalus_volans/mCynVol1/genomic_data/pacbio_hifi/m64330e_211212_152109.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Cynocephalus_volans/mCynVol1/genomic_data/pacbio_hifi/m64330e_211214_003514.hifi_reads.fastq.gz
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201003_145620.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201004_152428.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201014_175437.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201015_180930.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi

