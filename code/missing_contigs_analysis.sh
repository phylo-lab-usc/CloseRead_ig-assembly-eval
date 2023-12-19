source /etc/profile.d/modules.sh
module load conda
module load gcc/11.3.0
module load samtools/1.17
conda init
source /spack/conda/miniconda3/4.12.0/etc/profile.d/conda.sh
conda activate /home1/zhuyixin/.conda/envs/assembly

cd /home1/zhuyixin/sc2/ImmAssm/aligned_bam

# Filenames for the two BAM files
bam_file_1="/home1/zhuyixin/sc2/ImmAssm/aligned_bam/mCanLor1.missing.pri.bam"
bam_file_2="/home1/zhuyixin/sc1/ImmAssm/aligned_bam/combined/mCanLor1/mCanLor1_merged_sorted.pri.bam"
output_file_name="/home1/zhuyixin/sc2/ImmAssm/aligned_bam/mCanLor1.missing.matching_reads.sam"

temp_file="read_names.txt"

samtools view "$bam_file_1" | cut -f1 > "$temp_file"

# Use awk to filter reads from the second BAM file and print progress
samtools view -@ 60 "$bam_file_2" | awk -v RS='\n' -v FS='\t' -v OFS='\t' '
    NR==FNR {reads[$1]; next}
    FNR%1000000 == 0 {print "Processed " FNR " lines from the second file" > "/dev/stderr"}
    ($1 in reads)
' "$temp_file" - > "$output_file_name"


