#!/bin/bash

HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality
file_path="${HOME}/mammalianIG.txt"
while getopts 's:h:p:a:' option
do
    case "$option" in
        s) species=${OPTARG};;
        h) haploid=${OPTARG};;
        p) primary_result=${OPTARG};;
        a) alternate_result=${OPTARG};;
    esac
done

IggenePos_IG="${HOME}/gene_position/${species}/Ig_loci_details.txt"
output_file="${HOME}/gene_position/${species}/ref_loci_details.txt"

touch $IggenePos_IG
genes=("IGH" "IGK" "IGL")
if [ $alternate_result == "false" ] && [ $haploid == "False" ]
then
    while [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/${species}.alt.txt" ]; do
        #echo "Waiting for file $FILE_PATH to appear..."
        sleep 120  # Wait for 120 seconds before checking again
    done
    code/geneLociAutomated.sh -s ${species} -g alt
    for gene in "${genes[@]}"; do
        bed="gene_position/${species}/alt/${species}_${gene}_pos.sorted.bed"
        awk -v s="$species" -v p="alternate" -v g="$gene" '{
            # Increment the count for this chromosome
            count[$1]++
            # Update the start and end positions for this chromosome
            if (!start[$1] || $2 < start[$1]) {
                start[$1] = $2
            }
            if (!end[$1] || $3 > end[$1]) {
                end[$1] = $3
            }
        } END {
            # Iterate over the chromosomes and print the ones with at least 2 appearances
            for (chr in count) {
                if (count[chr] >= 5) {
                    print s, p, g, chr, start[chr], end[chr]
                }
            }
        }' "$bed" >> ${IggenePos_IG}
    done
fi
if [ $primary_result == "false" ]
then
    while [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/${species}.pri.txt" ]; do
        #echo "Waiting for file $FILE_PATH to appear..."
        sleep 120  # Wait for 120 seconds before checking again
    done
    code/geneLociAutomated.sh -s ${species} -g pri
    for gene in "${genes[@]}"; do
        bed="gene_position/${species}/pri/${species}_${gene}_pos.sorted.bed"
        awk -v s="$species" -v p="primary" -v g="$gene" '{
            # Increment the count for this chromosome
            count[$1]++
            # Update the start and end positions for this chromosome
            if (!start[$1] || $2 < start[$1]) {
                start[$1] = $2
            }
            if (!end[$1] || $3 > end[$1]) {
                end[$1] = $3
            }
        } END {
            # Iterate over the chromosomes and print the ones with at least 2 appearances
            for (chr in count) {
                if (count[chr] >= 5) {
                    print s, p, g, chr, start[chr], end[chr]
                }
            }
        }' "$bed" >> ${IggenePos_IG}
    done 
fi
touch ${IggenePos_IG}
echo "Igdetective finished"

final="${HOME}/gene_position/${species}/final.Ig_loci.txt"
rm -rf $final
touch "$final"
if [ $haploid == "False" ]
then
    haplotypes=("primary" "alternate")
else
    haplotypes=("primary")
fi
for gene in "${genes[@]}"; do
    # Iterate through each haplotype
    for hap in "${haplotypes[@]}"; do
        echo $hap
        # Check if the haplotype has an entry in refout
        found_in_ref=$(awk -v gene="$gene" -v hap="$hap" '$2 == hap && $3 == gene' "$output_file")
        if [[ ! -z "$found_in_ref" ]]; then
            # If found in refout, output those entries
            echo "$found_in_ref" >> "$final"
        else
            # If not found in refout, find and output from Igout
            echo $gene $hap
            echo $(awk -v gene="$gene" -v hap="$hap" '$2 == hap && $3 == gene' "${IggenePos_IG}")
            awk -v gene="$gene" -v hap="$hap" '$2 == hap && $3 == gene' "${IggenePos_IG}" >> "$final"
        fi
    done
done
echo "all finished"
