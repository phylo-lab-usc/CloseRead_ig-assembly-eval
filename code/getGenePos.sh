#!/bin/bash

#Count each gene's length
awk -F$'\t' '{print $5}' ../../mammalian_igdetective_v2.0/mCanLor1_igdetective/combined_genes_IGH.txt | awk -v FS="" '{print NF;}' > mCanLor1genes_IGH_length.txt 
#add a header "Length"
sed -i  "1s/.*/Length/" mCanLor1genes_IGH_length.txt 
#get Chromosome, start, strand infomation
awk -F$'\t' '{print $2 "\t" $3 "\t" $4}' combined_genes_IGH.txt > gene_IGH_start.txt
#merge length with the chromosome information
paste -d"\t" mCanLor1genes_IGH_length.txt gene_IGH_start.txt > temp1.txt 
#add start and length to get gene end position
awk -F"\t" 'NR==1 {val = "End"} NR > 1 {val = $1+$3} {$(++NF) = val; print}' OFS='\t' temp1.txt > temp2.txt
#remove the Length column
cut -f2- temp2.txt > gene_IGH_pos.txt 
#remove temporary files
rm -rf temp*.txt

## Next we want to convert gene_IGH_pos.txt into a bed file
# remove the headers
sed -i '1d' gene_IGH_pos.txt
# re-order file to convert into bed format
awk -F'\t' '{print $1 "\t" $2 "\t" $4 "\t" "." "\t" "0" "\t" $3 }' gene_IGH_pos.txt > gene_IGH_pos.bed

## Extract reads from BAM file
samtools view -b -L gene_position/mCanLor1/gene_IGH_pos.bed aligned_bam/mCanLor1/m64016_200910_161534_mapping.bam > extracted_bam/mCanLor1/m64016_200910_161534_extracted.bam