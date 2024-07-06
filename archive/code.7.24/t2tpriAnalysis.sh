#!/bin/bash
set -e -x

reffn=/home/egenge01/projects/IGL_ref_mod/reference_ready/modified_reference_renamed.fasta
IG_loci=/home/egenge01/projects/12_sample_test/IG_loci.bed
masked_ref=/home/egenge01/projects/12_sample_test/reference_IGloci_masked.fasta
scratch=$PWD
mask_ref=${scratch}/ref_IG_masked.fasta

function get_read_support_vdj3 {
    while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import ccs_bam
    do
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        bam_file="${scratch}/read_support/${sample}/ccs_to_pers/output.sorted.bam" # ccs reads to personalized reference
        ref="${scratch}/read_support/${sample}/ccs_to_pers/pers_ref.fasta" # personalized reference

        if [ ! -f "${bam_file}.bai" ]; then
            samtools index "$bam_file"
        fi

        for gene_type in "chr2" "chr22" "igh" #"ighc"
        do
            import_out="${base_outd}/${gene_type}/${sample}_make_gene_file_imported.csv"

            if [[ -f "$import_out" ]]; then
		modified_import_out="$import_out"
		
                tmp_file="${import_out}_read_support.tmp"
                echo "Total_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Percent_Accuracy,Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match" > "$tmp_file"
                header=$(head -n 1 "$import_out")
                IFS=',' read -ra header_cols <<< "$header"
                for i in "${!header_cols[@]}"; do
                    case "${header_cols[$i]}" in
                        "contig") contig_col=$i ;;
                        "REGION_start") start_col=$i ;;
                        "REGION_end") end_col=$i ;;
			"gene") gene_col=$i ;;  # Added case for 'gene'
                    esac
                done
		tmp_counts="${tmp_file}_counts"  # Temporary file to hold counts
		> "$tmp_counts"  # Clear or create the temp file for counts

                tail -n +2 "$import_out" | while IFS=, read -ra line
                do
                    contig="${line[$contig_col]}"
                    start=$(echo "${line[$start_col]}" | awk '{printf "%.0f", $1}')
                    end=$(echo "${line[$end_col]}" | awk '{printf "%.0f", $1}')
                    
		    gene="${line[$gene_col]}"  # Extract the 'gene' value using the identified column index
                    region="${contig}:${start}-${end}"

                    contig_filename=$(echo "$contig" | tr '/' '_')
                    tmp_bam="${base_outd}/${gene_type}/${contig_filename}_${start}_${end}.bam"
                    mkdir -p "$(dirname "$tmp_bam")"
                    samtools view -F 0x100 -F 0x800 -b "$bam_file" -o "$tmp_bam" -U "/dev/null" "${contig}:${start}-${end}"
                    samtools index "$tmp_bam"

		    samtools mpileup -f "$ref" -r "$region" "$tmp_bam" | \
			awk -v total_positions="$((end - start + 1))" -v sample="$sample" \
			'BEGIN {
    total_reads=0; mismatched_positions=0; matched_positions=0; positions_with_10x=0;
    mismatch_list=""; match_list="";
}
{
    total_reads += length($5);
    mismatches = length(gensub(/[.,]/, "", "g", $5));
    matches = length(gensub(/[^.,]/, "", "g", $5));
    mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches);
    match_list = (match_list == "" ? matches : match_list ":" matches);

    coverage = length($5);
    if (coverage >= 10) {
        positions_with_10x++;
    }

    mismatch_rate = mismatches / coverage;
    match_rate = matches / coverage;

    if (mismatch_rate > 0.2) {
        mismatched_positions++;
    }

    if (match_rate > 0.8) {
        matched_positions++;
    }
}
END {
    avg_reads_per_position = (total_positions > 0) ? total_reads / total_positions : 0;
    percent_accuracy = (matched_positions / total_positions) * 100;
    print total_positions, avg_reads_per_position, mismatched_positions, matched_positions, mismatch_list, match_list, percent_accuracy, positions_with_10x;}' OFS=',' >> "${tmp_file}_awk_out"
		    python match_subsequences3.py "$tmp_bam" "$contig" "$start" "$end" "$gene" "$import_out" > "${tmp_file}_py_out"
		    wait
		    paste -d ',' "${tmp_file}_awk_out" "${tmp_file}_py_out" >> "$tmp_file"
		    rm "${tmp_file}_awk_out" "${tmp_file}_py_out"
                    rm "$tmp_bam" "${tmp_bam}.bai"
                done

# New block to merge data and update Subject and Sample_Name
		if [[ -f "$import_out" && -f "$tmp_file" ]]; then
		    combined_file="${import_out%.csv}_combined.csv"
		    final_output="${import_out%.csv}_with_read_support.csv"
		    
		    # Merge the original import file with the new tmp file
		    paste -d ',' "$import_out" "$tmp_file" > "$combined_file"
		    rm -f "$tmp_file"
		    
		    # Update the Subject and Sample_Name columns in the final output
		    awk -v sample="$sample" 'BEGIN{FS=OFS=","} {
    if (NR == 1) {
        print;
    } else {
        $2 = sample;  # Update the 2nd column with $sample
        $3 = sample;  # Update the 3rd column with $sample
        print;
    }
}' "$combined_file" > "$final_output"

		    rm -f "$combined_file"
		fi
	    fi
        done
    done < extended_samples_paths.fofn
}