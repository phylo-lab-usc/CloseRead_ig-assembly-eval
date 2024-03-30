#!/bin/bash

HOME=/home1/zhuyixin/zhuyixin_proj/AssmQuality
file_path="${HOME}/mammalianIG.txt"
while getopts 's:h:p:a:' option
do
    case "$option" in
        s) species=${OPTARG};;
        h) haploid=${OPTARG};;
        p) pri_outdir=${OPTARG};;
        a) alt_ourdir=${OPTARG};;
    esac
done

check_and_print_loci() {
    local subject=$1
    local haplotype_prefix=$2
    local found_all_loci=true
    for locus in IGH IGK IGL; do
        # Find entries in the file that match the subject and locus
        read contig startPos endPos <<< $(grep -P "^${subject}\t" "$file_path" | grep -P "\t${locus}\t" | awk -v haplotype="$haplotype_prefix" '{ if ($2 == haplotype) print $4, $5, $6 }')
        if [[ -z $contig ]]; then
            found_all_loci=false
            break
        fi
    done
    echo "$found_all_loci"
}

save_loci_details() {
    local subject=$1
    local haplotype_prefix=$2
    local line_output=""
    for locus in IGH IGK IGL; do
        local grep_output=$(grep -P "^${subject}\t" "$file_path")
        line_output=$(echo "$grep_output" | grep -P "\t${locus}\t" | \
        awk -v haplotype="$haplotype_prefix" -v subject="$subject" '{
            if ($2 == haplotype) { 
                contig=$4; 
                startPos=$5; 
                endPos=$6; 
                if (haplotype == "alternate") contig="alt_"$4; 
                printf "%s\t%s\t%s\t%s\t%s\t%s\n", subject, haplotype, $3, contig, startPos, endPos;
            } }')
        if [[ -n $line_output ]]; then
            echo "$line_output" >> $output_file
        fi
    done
}

primary_result=$(check_and_print_loci "$species" "primary")
alternate_result=$(check_and_print_loci "$species" "alternate")
output_file="${HOME}/gene_position/${species}/ref_loci_details.txt"
prigenome="${HOME}/assemblies/${species}.pri.fasta"
altgenome="${HOME}/assemblies/${species}.alt.fasta"
IggenePos_IG="${HOME}/gene_position/${species}/Ig_loci_details.txt"

rm -rf ${output_file}
if [ $haploid == "False" ]
then
    save_loci_details "${species}" "alternate"
fi
save_loci_details "${species}" "primary"
echo "ref finished"

rm -rf ${IggenePos_IG}
if [ $alternate_result == "false" ] && [ $haploid == "False" ]
then
    sbatch --partition=gpu code/igDetective.sh ${altgenome} $alt_ourdir ${species} alt
elif [ $primary_result == "false"]
then
    sbatch --partition=gpu code/igDetective.sh ${prigenome} $pri_outdir ${species} pri
fi 
genes=("IGH" "IGK" "IGL")
if [ $alternate_result == "false" ] && [ $haploid == "False" ]
then
    while [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/${species}.alt.txt" ]; do
        echo "Waiting for file $FILE_PATH to appear..."
        sleep 30  # Wait for 30 seconds before checking again
    done
    code/geneLociAutomated.sh -s ${species} -g alt
    for gene in "${genes[@]}"; do
        bed="gene_position/${species}/alt/${species}_${gene}_pos.sorted.bed"
        awk -v s="$species" -v p="altrnate" -v g="$gene" '{
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
                if (count[chr] >= 2) {
                    print s, p, g, chr, start[chr], end[chr]
                }
            }
        }' "$bed" >> ${IggenePos_IG}
    done
elif [ $primary_result == "false"]
then
    while [ ! -f "/home1/zhuyixin/zhuyixin_proj/AssmQuality/igGene/${species}.pri.txt"]; do
        echo "Waiting for file $FILE_PATH to appear..."
        sleep 30  # Wait for 30 seconds before checking again
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
                if (count[chr] >= 2) {
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
if [ $haploid == "False" ]
then
    haplotypes=("primary" "alternate")
else
    haplotypes=("primary")
fi
for gene in "${genes[@]}"; do
    # Iterate through each haplotype
    for hap in "${haplotypes[@]}"; do
        # Check if the haplotype has an entry in refout
        found_in_ref=$(awk -v gene="$gene" -v hap="$hap" '$2 == hap && $3 == gene' "$output_file")
        if [ ! -z "$found_in_ref" ]; then
            # If found in refout, output those entries
            echo "$found_in_ref" >> "$final"
        else
            # If not found in refout, find and output from Igout
            awk -v gene="$gene" -v hap="$hap" '$2 == hap && $3 == gene' "${IggenePos_IG}" >> "$final"
        fi
    done
done
touch $final
