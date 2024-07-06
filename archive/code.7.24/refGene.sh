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
output_file="${HOME}/gene_position/${species}/ref_loci_details.txt"
touch $output_file

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
            } 
        }')
        if [[ -n $line_output ]]; then
            echo "$line_output" >> "$output_file"
        fi
    done
}

primary_result=$(check_and_print_loci "$species" "primary")
alternate_result=$(check_and_print_loci "$species" "alternate")

prigenome="${HOME}/assemblies/${species}.pri.fasta"
altgenome="${HOME}/assemblies/${species}.alt.fasta"
IggenePos_IG="${HOME}/gene_position/${species}/Ig_loci_details.txt"

if [ $haploid == "False" ]
then
    save_loci_details "${species}" "alternate"
fi
save_loci_details "${species}" "primary"

echo "$primary_result"
echo "$alternate_result"
rm -rf ${IggenePos_IG}

