import os
import subprocess

def extract_fasta_region(fasta_file, region, output_file):
    """Use samtools to extract a specific region from a fasta file."""
    cmd = f"samtools faidx {fasta_file} {region} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

# Path to the lociMeta.txt and assemblies directory
loci_meta_path = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/lociMeta.txt"
assemblies_dir = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/"
output_dir = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/IG_new"  # Output directory inside assemblies

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Read the lociMeta.txt file
with open(loci_meta_path, 'r') as loci_file:
    for line in loci_file:
        parts = line.strip().split()
        species = parts[0]
        assembly_type = parts[1].lower()  # 'primary' or 'alternate'
        ig_type = parts[2]  # 'IGH', 'IGK', 'IGL'
        chrom = parts[3]
        start = parts[4]
        end = parts[5]
        region = f"{chrom}:{start}-{end}"
        
        # Find the corresponding .merged.fasta file
        fasta_file = f"{species}.merged.fasta"
        fasta_path = os.path.join(assemblies_dir, fasta_file)
        
        # Output file path (distinguished by assembly type and placed in the IG directory)
        output_fasta = f"{species}_{assembly_type}_{ig_type}.fasta"
        output_path = os.path.join(output_dir, output_fasta)
        
        # Extract the region using samtools
        if os.path.exists(fasta_path):  # Ensure the fasta file exists
            extract_fasta_region(fasta_path, region, output_path)
        else:
            print(f"File not found: {fasta_path}")

print("Extraction complete.")
