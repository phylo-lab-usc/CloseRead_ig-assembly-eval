#!/bin/bash
#SBATCH --job-name=download    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

#module load python
#module spider anaconda
#module load conda
#conda init
#source /spack/apps/anaconda3/2021.05/etc/profile.d/conda.sh
#conda activate igortest

mkdir /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mMelMel3
cd /home1/zhuyixin/sc1/ImmAssm/hifi_fastq/mMelMel3

wget https://genomeark.s3.amazonaws.com/species/Meles_meles/mMelMel3/genomic_data/pacbio_hifi/m64094_200910_173211.ccs.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Meles_meles/mMelMel3/genomic_data/pacbio_hifi/m64094_200911_174739.ccs.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Meles_meles/mMelMel3/genomic_data/pacbio_hifi/m64125_200903_154937.ccs.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Meles_meles/mMelMel3/genomic_data/pacbio_hifi/m64125_201113_134826.ccs.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Meles_meles/mMelMel3/genomic_data/pacbio_hifi/m64125_201114_141552.ccs.bc1022_BAK8B_OA--bc1022_BAK8B_OA.bam
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201015_180930.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64016_200910_161534.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64089_200917_151241.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201003_145620.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201004_152428.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201014_175437.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi
#wget https://genomeark.s3.amazonaws.com/species/Canis_lupus_orion/mCanLor1/genomic_data/pacbio_hifi/m64094_201015_180930.ccs.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam.pbi

