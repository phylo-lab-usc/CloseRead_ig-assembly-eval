#contigs_to_keep = {"alt_scaffold_9", "scaffold_8"}
contigs_to_keep = {"contig_17636", "contig_74247", "contig_pri_8"}
genes_to_keep = {"V"}
# Define the position range
alt_start_position = 98910
alt_end_position = 1634400
pri_start_position = 75525176
pri_end_position = 77039516

with open('/home1/zhuyixin/sc1/AssmQuality/igGene/mCanLor1.lja.igdetective/combined_genes_IGH.txt', 'r') as input_file, open('/home1/zhuyixin/sc1/AssmQuality/igGene/mCanLor1.lja.igdetective/combined_genes_IGH.fasta', 'w') as output_file:
    # Skip the header line
    next(input_file)
    i = 0 
    for line in input_file:
        if not line.strip():  # Skip empty lines
            continue
        gene_type, contig, pos, strand, sequence, _, _ = line.strip().split('\t')
        if contig in contigs_to_keep:
            if gene_type in genes_to_keep:
                """
                #Convert pos to an integer for comparison
                pos = int(pos)
                # Determine the position range based on the contig name
                start_range = 0
                end_range = 0
                if "alt" in contig:
                    start_range = alt_start_position
                    end_range = alt_end_position
                else:
                    start_range = pri_start_position
                    end_range = pri_end_position
                invert = ""
                if start_range <= pos and pos <= end_range:
                    invert = "T"
                else:
                    invert = "F"
                # Create a FASTA entry with a header and sequence
                if "alt" in contig:
                    contig = "alt"
                else:
                    contig = "pri"
                header = f'>{i}_{gene_type}_{contig}{invert}_{pos}_{strand}'
                fasta_entry = f'{header}\n{sequence}'
                """
                if "contig_pri_8" in contig:
                    contig_label = "1pri"
                elif "contig_74247" in contig:
                    contig_label = "1Malt"
                else:
                    contig_label = "1alt"
                # Create the header for the FASTA enstry
                header = f'>{i}_{gene_type}_{contig_label}_{pos}_{strand}'
                
                # Create the FASTA entry
                fasta_entry = f'{header}\n{sequence}'
                i+=1
                # Write the FASTA entry to the output file
                output_file.write(fasta_entry + '\n')

