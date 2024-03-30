from Bio import SeqIO
import os

def consolidate_fasta(fasta_file, output_file, label_mapping_file):
    seq_records = {}
    contig_labels = {}

    # Parsing the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        desc = record.description.split('_')
        gene_type = desc[1]
        contig_label = desc[2]
        unique_id = gene_type + record.seq

        if unique_id in seq_records:
            contig_labels[unique_id].add(contig_label)
        else:
            seq_records[unique_id] = record
            contig_labels[unique_id] = {contig_label}

    # Writing unique sequences to a new FASTA file
    with open(output_file, 'w') as output_handle:
        for unique_id, record in seq_records.items():
            labels = ",".join(contig_labels[unique_id])
            record.id = record.id + "_" + labels
            record.description = ""  # Clearing the original description
            SeqIO.write(record, output_handle, "fasta")

    # Optionally, save the contig labels mapping
    with open(label_mapping_file, 'w') as mapping_handle:
        for unique_id, labels in contig_labels.items():
            mapping_handle.write(f"{unique_id}: {','.join(labels)}\n")

# Usage
consolidate_fasta("your_input_file.fasta", "unique_sequences.fasta", "contig_labels_mapping.txt")
