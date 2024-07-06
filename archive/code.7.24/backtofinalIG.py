import os

def read_and_write_gene_data(file_path, base_dir):
    """
    Read gene information from a file and write the classified information to species-specific files.
    """
    gene_info = {}

    # Open the file and process line by line
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            species, assembly, gene_type, chrom, start, end = parts

            # Initialize the data structure if species not encountered before
            if species not in gene_info:
                gene_info[species] = {
                    'IGH': [],
                    'IGK': [],
                    'IGL': []
                }
            
            # Format the region information
            region = f"{assembly} {gene_type} {chrom} {start} {end}"
            
            # Add the region to the correct list in the dictionary
            if gene_type in gene_info[species]:
                gene_info[species][gene_type].append(region)

    # Write the data to files based on species
    for species, types in gene_info.items():
        species_dir = os.path.join(base_dir, species)
        os.makedirs(species_dir, exist_ok=True)
        file_path = os.path.join(species_dir, 'final.Ig_loci.txt')
        with open(file_path, 'w') as out_file:
            for gene_type, regions in types.items():
                for region in regions:
                    out_file.write(f"{species} {region}\n")

def main():
    file_path = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/lociMeta.txt"  # Path to the file containing gene data
    base_dir = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/gene_position"  # Base directory for writing output
    
    # Process the gene data and write to respective files
    read_and_write_gene_data(file_path, base_dir)

if __name__ == "__main__":
    main()
