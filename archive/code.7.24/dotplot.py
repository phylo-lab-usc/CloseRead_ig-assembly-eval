import os
import subprocess
from collections import defaultdict

def run_gepard(x, y, species, gene_type, output_dir):
    """Run Gepard command line to generate dot plots."""
    os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists
    x_assembly_type = os.path.basename(x).split('_')[1]
    y_assembly_type = os.path.basename(y).split('_')[1] 
    output_file = f"{output_dir}/{species}.{gene_type}.{x_assembly_type}.{y_assembly_type}.png"
    cmd = f"java -cp /home1/zhuyixin/zhuyixin_proj/AssmQuality/code/Gepard-2.1.jar org.gepard.client.cmdline.CommandLine -seq {x} {y} -matrix /home1/zhuyixin/zhuyixin_proj/AssmQuality/code/edna.mat -outfile {output_file} -maxwidth 700 -maxheight 700"
    subprocess.run(cmd, shell=True, check=True)

def get_fasta_files(input_dir):
    """Organize fasta files by assembly type, species, and gene type."""
    fasta_files = {'primary': defaultdict(lambda: defaultdict(list)),
                   'alternate': defaultdict(lambda: defaultdict(list))}
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fasta'):
            parts = file_name.split('_')
            species = parts[0]
            assembly_type = parts[1]  # 'primary' or 'alternate'
            gene_type = parts[2].replace('.fasta', '')
            fasta_files[assembly_type][species][gene_type].append(os.path.join(input_dir, file_name))
    return fasta_files

def process_comparisons(fasta_files, output_dir_base):
    """Process each species and gene type to run Gepard comparisons."""
    for species in set(fasta_files['primary'].keys()).union(fasta_files['alternate'].keys()):
        for gene_type in set(fasta_files['primary'][species].keys()).union(fasta_files['alternate'][species].keys()):
            primary_files = fasta_files['primary'][species].get(gene_type, [])
            alternate_files = fasta_files['alternate'][species].get(gene_type, [])
            output_dir = os.path.join(output_dir_base, species)

            # Primary vs Primary
            for i in range(len(primary_files)):
                run_gepard(primary_files[i], primary_files[i], species, gene_type, output_dir)

            # Alternate vs Alternate
            for i in range(len(alternate_files)):
                run_gepard(alternate_files[i], alternate_files[i], species, gene_type, output_dir)

            # Primary vs Alternate
            for primary_file in primary_files:
                for alternate_file in alternate_files:
                    run_gepard(primary_file, alternate_file, species, gene_type, output_dir)

input_dir = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/IG_new"
output_dir_base = "/home1/zhuyixin/zhuyixin_proj/AssmQuality/errorPlots"
fasta_files = get_fasta_files(input_dir)
process_comparisons(fasta_files, output_dir_base)
print("Dot plot generation complete.")
