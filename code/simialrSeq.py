from Bio import AlignIO
from Bio import pairwise2
from Bio.Seq import Seq
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
import numpy as np

def calculate_similarity(seq1, seq2):
    # Ensure sequences are of the same length for comparison
    assert len(seq1) == len(seq2), "Sequences must be of the same length."

    # Count mismatches
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    
    # Calculate similarity
    similarity = 1 - (mismatches / len(seq1))
    return similarity

def create_similarity_matrix(alignment):
    num_sequences = len(alignment)
    similarity_matrix = np.zeros((num_sequences, num_sequences))
    for i in range(num_sequences):
        print(i)
        for j in range(i+1, num_sequences):
            similarity = calculate_similarity(alignment[i].seq, alignment[j].seq)
            similarity_matrix[i][j] = similarity_matrix[j][i] = similarity
    return similarity_matrix

def cluster_sequences(similarity_matrix, threshold=0.0001):
    # Convert similarity matrix to distance matrix
    distance_matrix = 1 - similarity_matrix
    np.fill_diagonal(distance_matrix, 0)

    # Ensure distance matrix is condensed, keeping only the upper triangle
    condensed_distance_matrix = squareform(distance_matrix)
    # Perform hierarchical clustering
    Z = linkage(condensed_distance_matrix, 'average')
    # Form flat clusters
    clusters = fcluster(Z, t=threshold, criterion='distance')
    return clusters

# Example usage
file_path = '/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/combined_genes_IGH.phy'
alignment = AlignIO.read(file_path, "phylip")
print("read complete")
similarity_matrix = create_similarity_matrix(alignment)
print(similarity_matrix)
print("similarity_matrix complete")
clusters = cluster_sequences(similarity_matrix)
print("cluster complete")


# Assuming 'clusters' is your array of cluster IDs from fcluster
# And 'alignment' is the BioPython MultipleSeqAlignment object

for i, cluster_id in enumerate(clusters):
    sequence_id = alignment[i].id  # Access the ID of the i-th sequence in the alignment
    print(f"{sequence_id} belongs to cluster {cluster_id}")

import csv

# Assuming 'clusters' contains your cluster IDs and 'alignment' is your sequence alignment
with open('/home1/zhuyixin/sc1/AssmQuality/treeBuilding/mCanLor/IGH.cluster_assignments.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Sequence ID', 'Cluster ID'])  # Write the header
    
    for i, cluster_id in enumerate(clusters):
        sequence_id = alignment[i].id  # Get the sequence ID
        writer.writerow([sequence_id, cluster_id])
