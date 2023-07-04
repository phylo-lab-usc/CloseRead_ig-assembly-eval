#!/bin/bash

awk -F$'\t' '{print $5}' ../../mammalian_igdetective_v2.0/mCanLor1_igdetective/combined_genes_IGH.txt | awk -v FS="" '{print NF;}' > mCanLor1genes_IGH_length.txt 
sed -i  "1s/.*/Length/" mCanLor1genes_IGH_length.txt 
awk -F$'\t' '{print $2 "\t" $3 "\t" $4}' combined_genes_IGH.txt > gene_IGH_start.txt
paste -d"\t" mCanLor1genes_IGH_length.txt gene_IGH_start.txt > temp1.txt 
awk -F"\t" 'NR==1 {val = "End"} NR > 1 {val = $1+$3} {$(++NF) = val; print}' OFS='\t' temp1.txt > temp2.txt
cut -f2- temp2.txt > gene_IGH_pos.txt 
rm -rf temp*.txt