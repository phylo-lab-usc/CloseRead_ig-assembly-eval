import os
from collections import defaultdict
import argparse

# Function to get count of V, D, J genes
def count_gene_types(file_path):
    counts = {'V': 0, 'D': 0, 'J': 0}
    with open(file_path, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            if line.strip():  # Check if line is not empty
                gene_type = line.split()[0]
                if gene_type in counts:
                    counts[gene_type] += 1
    return counts

def count_rows(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file if line.strip())
    
def process_chromosomes(file_path):
    """
    Processes BED file and organizes it into a dictionary.

    Parameters:
    file_path (str): The path to the BED file.

    Returns:
    dict: A dictionary where keys are chromosome identifiers (e.g., 'chr1', 'chr2') and values are lists of positions.
    """
    chrom_data = defaultdict(list)
    with open(file_path, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            parts = line.strip().split()
            chrom = parts[1]
            start = int(parts[2])
            chrom_data[chrom].append(start)
    return chrom_data

def analyze_chromosome_occurrences(chrom_data):
    """
    Analyzes chromosome occurrences to determine which chromosome has the most positions and the range of those positions.

    Parameters:
    chrom_data (dict): A dictionary where keys are chromosome identifiers (e.g., 'chr1', 'chr2') and values are lists of positions.

    Returns:
    tuple: A tuple containing the chromosome with the most occurrences, the number of occurrences, the minimum position, and the maximum position.
           If no chromosomes are present, returns (None, 0, None, None).
    """
    most_occurrences = 0
    target_chrom = None
    for chrom, positions in chrom_data.items():
        if len(positions) > most_occurrences:
            most_occurrences = len(positions)
            target_chrom = chrom
    if target_chrom:
        return (target_chrom, most_occurrences, min(chrom_data[target_chrom]), max(chrom_data[target_chrom]))
    return (None, 0, None, None)

# Main processing function
def process_ig_genes(species, home):
    gene_types = ['IGH', 'IGK', 'IGL']
    assembly_types = {'pri': 'primary', 'alt': 'alternate'}
    base_dir = f'{home}/igGene'
    output_dir = f'{home}/gene_position'
    primary_data = None
    outfile = f"{output_dir}/{species}.final.Ig_loci.txt"
    if os.path.exists(outfile):
        os.remove(outfile)
        print(f"File {outfile} has been overwrote.")
    for gene in gene_types:
        primary_file_path = f"{base_dir}/{species}.pri.igdetective/combined_genes_{gene}.txt"
        # Check if the primary file exists and has more than one row
        if os.path.exists(primary_file_path) and count_rows(primary_file_path) > 1:
            primary_chrom_data = process_chromosomes(primary_file_path)
            # Analyze the chromosome data to find the chromosome with the most occurrences
            primary_data = analyze_chromosome_occurrences(primary_chrom_data)
        # Iterate over each assembly type (primary and alternate)
        for assembly_code, assembly_name in assembly_types.items():
            file_path = f"{base_dir}/{species}.{assembly_code}.igdetective/combined_genes_{gene}.txt"
            if os.path.exists(file_path) and count_rows(file_path) > 1:
                chrom_data = process_chromosomes(file_path)
                chrom_info = analyze_chromosome_occurrences(chrom_data)
                #process primary haplotype normally, while only process alternate if its locus length is >= 1/4 of primary locus length
                if chrom_info[0] and (assembly_code == 'pri' or (primary_data and chrom_info[3] - chrom_info[2] >= (primary_data[3] - primary_data[2])/ 4)):
                    #only write the chromosomes if it has occurance > 2 time (fiter out outliers)
                    if chrom_info[1] > 2:
                        with open(f"{output_dir}/{species}.final.Ig_loci.txt", 'a') as f:
                            f.write(f"{species} {assembly_name} {gene} {chrom_info[0]} {chrom_info[2]} {chrom_info[3]}\n")
                        # with open(f'{home}/haplotype.txt', 'a') as hfile:
                        #     hfile.write(f"{species} {assembly_name} {gene} True\n")
                # else:
                #     with open(f'{home}/haplotype.txt', 'a') as hfile:
                #         hfile.write(f"{species} {assembly_name} {gene} False\n")
                # counting IGH primary haplotype VDJ genes
                # if gene == 'IGH' and assembly_code == 'pri':
                #     counts = count_gene_types(file_path)
                #     with open(f'{home}/VgeneCount.txt', 'a') as vfile:
                #         vfile.write(f"{species} {assembly_name} {counts['V']}\n")
                #     with open(f'{home}/DgeneCount.txt', 'a') as dfile:
                #         dfile.write(f"{species} {assembly_name} {counts['D']}\n")
                #     with open(f'{home}/JgeneCount.txt', 'a') as jfile:
                #         jfile.write(f"{species} {assembly_name} {counts['J']}\n")
            # else:
            #     with open(f'{home}/haplotype.txt', 'a') as hfile:
            #         hfile.write(f"{species} {assembly_name} {gene} False\n")

# Usage
def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='final IG txt..')

    # Required arguments
    parser.add_argument('species', help='species name')
    parser.add_argument('home', help='home dir')
    args = parser.parse_args()

    process_ig_genes(args.species, args.home)



if __name__ == "__main__":
    main()