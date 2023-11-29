samtools view -@60 -b human_merged_sorted.bam NC_000014.9 > human.chr14.bam
samtools sort -@ 8 -o human.chr14.sam human.chr14.bam
samtools view -F 0x800 -F 0x100 -@ 60 human.chr14.bam > human.chr14.pri.sam

cat human.chr14.pri.sam | while read readID readFlag chr pos mapq cigar rnext pnext tlen seq qual tagMS tag
do
    long_blocks_sum=0
    long_blocks_sum=$(echo $cigar | grep -oE '[0-9]+M' | awk '{ split($0, a, "M"); if(a[1] >= 500) sum += a[1] } END { print sum+0 }')
    total_length=$(echo $cigar | grep -oE '[0-9]+[MDI]' | awk '{ split($0, a, /[MDI]/); sum += a[1] } END { print sum }')
    fraction=0
    if [ $total_length -ne 0 ]; then
        fraction=$(awk "BEGIN {print $long_blocks_sum/$total_length}")
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${readID}" "${chr}" "${pos}" "${mapq}" "${total_length}" "${fraction}" >> human.chr14.pri.long_blocks.txt
done



