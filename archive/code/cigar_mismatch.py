#!/usr/bin/env python3

import pysam
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import seaborn as sns

def calculate_mismatches(read):
    """ Calculate the number of mismatches using NM tag and CIGAR string. """
    nm_tag = read.get_tag('NM') if read.has_tag('NM') else 0

    insertions = sum(length for op, length in read.cigartuples if op == 1)  # Sum of 'I'
    deletions = sum(length for op, length in read.cigartuples if op == 2)  # Sum of 'D'
    ambiguous = sum(length for op, length in read.cigartuples if op == 3)  # Sum of 'N'

    mismatches = nm_tag - insertions - deletions - ambiguous
    return mismatches

def process_bam_file(bam_file_path, output_file_path):
    """ Process a BAM file to estimate mismatches for each read. """
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    with open(output_file_path, 'w') as file:
        for read in bamfile.fetch():
            if not read.is_unmapped:
                mismatches = calculate_mismatches(read)
                read_length = read.query_length
                mismatch_rate = mismatches / read_length
                read_name = read.query_name
                chromosome = bamfile.get_reference_name(read.reference_id)
                position = read.reference_start
                mapping_quality = read.mapping_quality
                file.write(f"{read_name}\t{chromosome}\t{position}\t{mapping_quality}\t{read_length}\t{mismatches}\t{mismatch_rate}\n")
    bamfile.close()

def process_sam_file(sam_file_path, output_file_path):
    """ Process a SAM file to estimate mismatches for each read. """
    samfile = pysam.AlignmentFile(sam_file_path, "r")  # Open in "r" mode for SAM file
    with open(output_file_path, 'w') as file:
        for read in samfile.fetch():
            if not read.is_unmapped:
                mismatches = calculate_mismatches(read)
                read_length = read.query_length
                mismatch_rate = mismatches / read_length
                read_name = read.query_name
                chromosome = samfile.get_reference_name(read.reference_id)
                position = read.reference_start
                mapping_quality = read.mapping_quality
                file.write(f"{read_name}\t{chromosome}\t{position}\t{mapping_quality}\t{read_length}\t{mismatches}\t{mismatch_rate}\n")
    samfile.close()

def plot_histograms(file_path, output_file):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
        sns.set(style="whitegrid")
        colors = ["#87CEEB", "#FDA50F", "#228B22", "#708090"]
        mismatch = np.loadtxt(file_path, usecols=[5])
        mismatch_rate = np.loadtxt(file_path, usecols=[6])

        # Create the histogram
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
        sns.histplot(mismatch, bins=10, kde=False, ax=axes[0], color=colors[0], edgecolor='black')
        sns.histplot(mismatch_rate, bins=10, kde=False, ax=axes[1], color=colors[1], edgecolor='black')

        # Add titles and labels with bold font
        axes[0].set_title('Histogram of Mismatches', fontsize=18, fontweight='bold')
        axes[0].set_xlabel('Number of Mismatches', fontsize=14, fontweight='bold')
        axes[0].set_ylabel('Frequency', fontsize=14, fontweight='bold')
        axes[0].tick_params(axis='both', which='major', labelsize=12)

        axes[1].set_title('Histogram of Mismatch Rate', fontsize=18, fontweight='bold')
        axes[1].set_xlabel('Mismatches / Read Length', fontsize=14, fontweight='bold')
        axes[1].set_ylabel('Frequency', fontsize=14, fontweight='bold')
        axes[1].tick_params(axis='both', which='major', labelsize=12)
     
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(output_file, format='png', dpi=300)

def count_mismatches(file_path):
    # Read the data
    mismatch_rate = np.loadtxt(file_path, usecols=[6])

    # Count reads based on mismatch criteria
    more_than_001 = sum(mismatch_rate > 0.01)
    less_or_equal_001 = sum(mismatch_rate <= 0.01)

    return more_than_001, less_or_equal_001

def plot_mismatch_bar(more_than_001, less_or_equal_001, output_file):
    # Data for plotting
    categories = ['> 0.01 Mismatches Rate', '≤ 0.01 Mismatches Rate']
    counts = [more_than_001, less_or_equal_001]

    # Create the bar plot
    sns.set_style("whitegrid")
    colors = ["#87CEEB", "#FDA50F", "#228B22", "#708090"]

    plt.figure(figsize=(8, 6))
    sns.barplot(x=categories, y=counts, palette=colors)
    plt.title('Mismatches Rate', fontsize=18, fontweight='bold')
    plt.xlabel('Mismatches Rate', fontsize=14, fontweight='bold')
    plt.ylabel('Frequency', fontsize=14, fontweight='bold')

    # Set bold font for the tick labels
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_file, format='png', dpi=300)

def main():
    if len(sys.argv) != 5:
        print("Usage: python cigar_mismatch.py <input_file_path> <output_file_path> <output_plot_file> <output_barplot_file>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    output_plot_file = sys.argv[3]
    output_barplot_file = sys.argv[4]
    if input_file_path.endswith('.bam'):
        process_bam_file(input_file_path, output_file_path)
    elif input_file_path.endswith('.sam'):
        process_sam_file(input_file_path, output_file_path)
    plot_histograms(output_file_path, output_plot_file)
    more_than_001, less_or_equal_001 = count_mismatches(output_file_path)
    print(f"More than 0.01: {more_than_001}")
    plot_mismatch_bar(more_than_001, less_or_equal_001, output_barplot_file)

if __name__ == "__main__":
    main()
