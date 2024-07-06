def extract_and_modify_sequence_names_with_header(phylip_file_path, output_file_path):
    with open(phylip_file_path, 'r') as file, open(output_file_path, 'w') as output_file:
        # Write the header line
        output_file.write("name\tbranch_color\tleaf_dot_color\tbar1_height\tbar1_color\tbar2_height\tbar2_color\n")
        
        header = file.readline()  # Read the header line from the PHYLIP file
        num_sequences, _ = map(int, header.split())  # Extract the number of sequences (ignore sequence length here)
        
        block_count = 0  # To identify the first block
        for line in file:
            if line.strip():  # Ignore empty lines
                if block_count < num_sequences:
                    # In the first block, extract names
                    sequence_name = line[:10].strip()  # Assuming names are within the first 10 characters
                    
                    # Determine the second column based on 'alt'
                    if '_alt' in sequence_name:
                        forth_column = '-1'
                        fifth_column = 'bp_orange'
                        sixth_column = '0'
                        seventh_column = 'bp_light_purple'
                    elif '_pri' in sequence_name:
                        forth_column = '-1'
                        fifth_column = 'bp_light_orange'
                        sixth_column = '0'
                        seventh_column = 'bp_light_purple'
                    elif '1Malt' in sequence_name:
                        forth_column = '1'
                        fifth_column = 'bp_green'
                        sixth_column = '1'
                        seventh_column = 'bp_light_purple'
                    elif '1alt' in sequence_name:
                        forth_column = '1'
                        fifth_column = 'bp_green'
                        sixth_column = '0'
                        seventh_column = 'bp_light_purple'
                    elif '1pri' in sequence_name:
                        forth_column = '1'
                        fifth_column = 'bp_light_green'
                        sixth_column = '0'
                        seventh_column = 'bp_light_purple'
                    
                    # Determine the third column based on specific strings
                    if 'T' in sequence_name:
                        second_column = 'bp_blue'
                        third_column = 'bp_blue'
                    elif 'F' in sequence_name:
                        second_column = 'bp_pink'
                        third_column = 'bp_pink'
                    else:
                        second_column = ''
                        third_column = ''  # Or some default value if needed
                    
                    # Write the line to the output file
                    output_file.write(f"{sequence_name}\t{third_column}\t{second_column}\t{forth_column}\t{fifth_column}\t{sixth_column}\t{seventh_column}\n")
                    block_count += 1
                # No need to process further once all names in the first block are collected
                if block_count == num_sequences:
                    break


extract_and_modify_sequence_names_with_header('/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/combined_genes_IGH.phy', '/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/combined_genes_IGH.txt')
