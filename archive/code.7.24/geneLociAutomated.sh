#!/bin/bash

HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality
while getopts 's:g:' option
do
    case "$option" in
        s) species=${OPTARG};;
        g) genome=${OPTARG};;
    esac
done

#create output directories
mkdir ${HOME}/gene_position/
mkdir ${HOME}/gene_position/${species}/
mkdir ${HOME}/gene_position/${species}/${genome}

for gene in IGH IGK IGL
do
    #Count each gene's length
    sed '1d' ${HOME}/igGene/${species}.${genome}.igdetective/combined_genes_${gene}.txt | awk -F$'\t' '{print $5}' | awk -v FS="" '{print NF;}' > ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_length.txt
    #get Chromosome, start, strand infomation
    sed '1d' ${HOME}/igGene/${species}.${genome}.igdetective/combined_genes_${gene}.txt | awk -F$'\t' '{print $2 "\t" $3 "\t" $4}' > ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_start.txt
    #expand start to get gene position
    awk -F"\t" 'NR >= 1 {val = $2-20000} {$(++NF) = val; print}' OFS='\t' ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_start.txt > ${HOME}/gene_position/${species}/${genome}/temp1.txt
    cut -f1,3,4 ${HOME}/gene_position/${species}/${genome}/temp1.txt > ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_pos.txt
    # re-order file to convert into bed format
    awk -F'\t' '{print $1 "\t" $3 "\t" $3+40000 "\t" NR "\t" "0" "\t" $2 }' ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_pos.txt > ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_pos.bed
    #sort the bed file
    sort -k1,1 -k2,2n -k3,3n ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_pos.bed > ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_pos.sorted.bed
    #replace negative corrdinates with 0
    sed -i 's/-[0-9][0-9]*/0/' ${HOME}/gene_position/${species}/${genome}/${species}_${gene}_pos.sorted.bed
    rm -rf ${HOME}/gene_position/${species}/${genome}/temp1.txt
done





    