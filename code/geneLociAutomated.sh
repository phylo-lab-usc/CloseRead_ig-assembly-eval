#!/bin/bash
#SBATCH --job-name=gene_automated    # Job name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                    # Run on a single CPU 
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=gene_automated.log   # Standard output and error log


#HOME=/home1/zhuyixin/sc1/ImmAssm
while getopts 'f:p:' option
do
    case "$option" in
        f) function=${OPTARG};;
        p) HOME=${OPTARG};;
    esac
done
mkdir ${HOME}/gene_position/${function}
cat ${HOME}/igv_name.txt | while read line
do 
    #create output directories
    mkdir ${HOME}/gene_position/${function}/${line}
    #Select only the productive gene
    awk -F$'\t' '$6 == "True"' ${HOME}/mammalian_igdetective_v2.0/${line}_igdetective/combined_genes_IGH.txt > ${HOME}/gene_position/${function}/${line}/${line}_genes_IGH_productive.txt
    #Select only the non-productive gene
    awk -F$'\t' '$6 == "False"' ${HOME}/mammalian_igdetective_v2.0/${line}_igdetective/combined_genes_IGH.txt > ${HOME}/gene_position/${function}/${line}/${line}_genes_IGH_nonproductive.txt
    #save as variable
    if [ "$function" = "functional" ]; then
        IGHloci=${line}_genes_IGH_productive.txt
    else
        IGHloci=${line}_genes_IGH_nonproductive.txt
    fi
    #Count each gene's length
    awk -F$'\t' '{print $5}' ${HOME}/gene_position/${function}/${line}/${IGHloci} | awk -v FS="" '{print NF;}' > ${HOME}/gene_position/${function}/${line}/${line}_genes_IGH_length.txt
    #add a header "Length"
    #sed -i  "1s/.*/Length/" ${HOME}/gene_position/${function}/${line}/${line}_genes_IGH_length.txt 
    #get Chromosome, start, strand infomation
    awk -F$'\t' '{print $2 "\t" $3 "\t" $4}' ${HOME}/gene_position/${function}/${line}/${IGHloci} > ${HOME}/gene_position/${function}/${line}/gene_IGH_start.txt
    #merge length with the chromosome information
    paste -d"\t" ${HOME}/gene_position/${function}/${line}/${line}_genes_IGH_length.txt ${HOME}/gene_position/${function}/${line}/gene_IGH_start.txt > ${HOME}/gene_position/${function}/${line}/temp1.txt 
    #add start and length to get gene end position
    awk -F"\t" 'NR >= 1 {val = $1+$3} {$(++NF) = val; print}' OFS='\t' ${HOME}/gene_position/${function}/${line}/temp1.txt > ${HOME}/gene_position/${function}/${line}/gene_IGH_pos.txt
    #remove the Length column
    #cut -f2- ${HOME}/gene_position/${function}/${line}/temp2.txt > ${HOME}/gene_position/${function}/${line}/gene_IGH_pos.txt 
    #remove temporary files
    rm -rf ${HOME}/gene_position/${function}/${line}/temp*.txt 
    ## Next we want to convert gene_IGH_pos.txt into a bed file
    # remove the headers
    #sed -i '1d' ${HOME}/gene_position/${function}/${line}/gene_IGH_pos.txt
    # re-order file to convert into bed format
    awk -F'\t' '{print $2 "\t" $3 "\t" $5 "\t" NR "\t" "0" "\t" $4 }' ${HOME}/gene_position/${function}/${line}/gene_IGH_pos.txt > ${HOME}/gene_position/${function}/${line}/gene_IGH_pos.bed
    #sort the bed file
    sort -k1,1 -k2,2n -k3,3n ${HOME}/gene_position/${function}/${line}/gene_IGH_pos.bed > ${HOME}/gene_position/${function}/${line}/gene_IGH_pos_sorted.bed
done

    