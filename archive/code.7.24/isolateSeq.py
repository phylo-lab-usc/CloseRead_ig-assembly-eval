from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Load the FASTA file
#fasta_sequences = SeqIO.to_dict(SeqIO.parse("/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/Zebra.fasta", "fasta"))
with open("/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/Zebra.fasta", "r") as fasta_file:
    single_sequence = next(SeqIO.parse(fasta_file, "fasta")).seq

extracted_sequences = []  # List to store extracted sequences

# Open your data file and extract the sequences
with open("/home1/zhuyixin/zhuyixin_proj/AssmQuality/zebra.igh.txt", "r") as file:
    next(file)
    for line in file:
        parts = line.strip().split()  # Assuming the file is tab-separated
        if len(parts) < 5:
            continue  # Skip lines that don't have enough columns

        locus, type, name, start, end = parts
        start, end = int(start), int(end)  # Convert to integers

        # Assuming 'locus' corresponds to the sequence ID in the FASTA file
        sequence = single_sequence[start-1:end]  # -1 for 0-based indexing

        # Create a SeqRecord object for the extracted sequence
        sequence_id = f"{locus}{type}-{name}"
        seq_record = SeqRecord(sequence, id=sequence_id, description="")
        extracted_sequences.append(seq_record)

# Save the extracted sequences to a new FASTA file
#SeqIO.write(extracted_sequences, "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/Zebra.igh.fasta", "fasta")
with open("/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/Zebra.igh.fasta", "w") as output_file:
    for seq_record in extracted_sequences:
        output_file.write(f">{seq_record.id}\n{str(seq_record.seq)}\n")
