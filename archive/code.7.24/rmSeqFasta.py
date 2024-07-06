from Bio import SeqIO

# Load sequences, excluding the one named "xxx"
sequences = SeqIO.parse("/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/mCynVol1.corr.merged.fasta", "fasta")
filtered_sequences = [seq for seq in sequences if seq.id != "atg001440l_1"]

# Write the remaining sequences to a new file
SeqIO.write(filtered_sequences, "/home1/zhuyixin/zhuyixin_proj/AssmQuality/assemblies/mCynVol1.corr.merged.fasta1", "fasta")
