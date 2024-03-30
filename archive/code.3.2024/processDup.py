import csv
from collections import defaultdict

def parse_fasta(fasta_file):
    """
    Parse the FASTA file to extract sequences, their associated metadata, and truncated names.
    Sequences are keyed by their full sequence for uniqueness, with metadata including a list
    of positions categorized by sample type and the first 10 characters of the original header for naming.
    """
    sequences = defaultdict(lambda: {"metadata": defaultdict(list), "name": ""})
    with open(fasta_file, 'r') as f:
        header_name = ""
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                header_name = header[:10]  # Truncate the ID part of the header to the first 10 characters
                _, _, sample_type,  position, _ = header.split('_')
                sample_category = categorize_sample(sample_type)
            else:
                seq = line.strip()
                sequences[seq]["name"] = header_name  # Assign truncated name
                if sample_category:
                    sequences[seq]["metadata"][sample_category].append(position)
    return sequences

def categorize_sample(sample_type):
    """
    Categorize the sample type into one of the four specified categories.
    """
    if '1alt' in sample_type or '1Malt' in sample_type:
        return '1alt'
    elif '1pri' in sample_type:
        return '1pri'
    elif 'altT' in sample_type or 'altF' in sample_type:
        return '2alt'
    elif 'priT' in sample_type or 'priF' in sample_type:
        return '2pri'
    return None

def write_outputs(sequences, output_fasta, output_metadata):
    """
    Write the unique sequences to a FASTA file using the first 10 characters of the header as the sequence name,
    and metadata for duplicated sequences to a CSV file with 'NA' for missing categories.
    """
    with open(output_fasta, 'w') as fasta_out, open(output_metadata, 'w', newline='') as meta_out:
        writer = csv.DictWriter(meta_out, fieldnames=['Sequence Name', '1alt', '1pri', '2alt', '2pri'])
        writer.writeheader()
        
        for seq, info in sequences.items():
            row = {key: ';'.join(value) if value else "NA" for key, value in info["metadata"].items()}
            row['Sequence Name'] = info["name"]
            
            if sum(1 for v in info["metadata"].values() if v) > 1:
                # Sequence is duplicated, write its metadata to CSV
                writer.writerow(row)
            else:
                # Writing sequence to FASTA with its truncated header name
                fasta_out.write(f">{info['name']}\n{seq}\n")

fasta_file = "/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/combined_genes_IGH.fasta"
output_fasta = "/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/IGH.unique_sequences.fasta"
output_metadata = "/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/IGH.duplicated_sequences_metadata.csv"

sequences = parse_fasta(fasta_file)
print(sequences)
write_outputs(sequences, output_fasta, output_metadata)
