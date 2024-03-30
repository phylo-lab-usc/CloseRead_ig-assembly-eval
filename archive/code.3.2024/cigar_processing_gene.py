#!/usr/bin/env python3

import pysam
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import seaborn as sns
from intervaltree import Interval, IntervalTree
import argparse
import os
import warnings


def load_bed_to_interval_tree(bed_file_path):
    # Dictionary to hold an interval tree for each chromosome
    trees = {}
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            trees[chrom].addi(start, end)
    return trees

def calculate_mismatches(read):
    """ Calculate the number of mismatches using NM tag and CIGAR string. """
    nm_tag = read.get_tag('NM') if read.has_tag('NM') else 0

    insertions = sum(length for op, length in read.cigartuples if op == 1)  # Sum of 'I'
    deletions = sum(length for op, length in read.cigartuples if op == 2)  # Sum of 'D'
    ambiguous = sum(length for op, length in read.cigartuples if op == 3)  # Sum of 'N'
    soft_clipping = sum(length for op, length in read.cigartuples if op == 4)  # Sum of 'S'
    hard_clipping = sum(length for op, length in read.cigartuples if op == 5)  # Sum of 'H'

    mismatches = nm_tag - insertions - deletions - ambiguous
    longindels = 0
    total_indel_length = 0
    for operation, length in read.cigartuples:
        # Insertion or Deletion longer than 2
        if (operation == 1 or operation == 2) and length > 2:
            longindels += 1
            total_indel_length += length
    return mismatches, longindels, total_indel_length, soft_clipping, hard_clipping

def process_bam_file(bam_file_path, igh_bed_file_path, igk_bed_file_path, igl_bed_file_path, output_dir):
    """ Process a BAM file to estimate mismatches for each read. """
    bed_trees = [
    load_bed_to_interval_tree(igh_bed_file_path),
    load_bed_to_interval_tree(igk_bed_file_path),
    load_bed_to_interval_tree(igl_bed_file_path)
    ]
    # Output file names
    output_file_names = [os.path.join(output_dir, "IGH.txt"), os.path.join(output_dir, "IGK.txt"), os.path.join(output_dir, "IGL.txt")]
    non_overlap_file_name = os.path.join(output_dir,"nonIG.txt")

    # Open output files
    output_files = [open(name, "w") for name in output_file_names]
    non_overlap_file = open(non_overlap_file_name, "w")

    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    for read in bamfile:
        if not read.is_unmapped:
            mismatches, longindels, total_indel_length, soft_clipping, hard_clipping = calculate_mismatches(read)
            read_length = read.query_length
            if mismatches != 0:
                mismatch_rate = mismatches / read_length
            else:
                mismatch_rate = 0
            read_name = read.query_name
            chromosome = bamfile.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            mapping_quality = read.mapping_quality
            indel_rate = total_indel_length / read_length
            # Check for overlap with each BED file's intervals
            found_overlap = False  # Flag to check if an overlap was found
            for i, tree_dict in enumerate(bed_trees):
                if chromosome in tree_dict and tree_dict[chromosome].overlaps(start, end):
                    output_files[i].write(f"{read_name}\t{chromosome}\t{start}\t{read_length}\t{mapping_quality}\t{mismatches}\t{mismatch_rate}\t{longindels}\t{total_indel_length}\t{indel_rate}\t{soft_clipping}\t{hard_clipping}\n")
                    found_overlap = True
                    break
            # If no overlap was found, write to the non-overlapping file
            if not found_overlap:
                non_overlap_file.write(f"{read_name}\t{chromosome}\t{start}\t{read_length}\t{mapping_quality}\t{mismatches}\t{mismatch_rate}\t{longindels}\t{total_indel_length}\t{indel_rate}\t{soft_clipping}\t{hard_clipping}\n")
    bamfile.close()
    for file in output_files:
        file.close()
    non_overlap_file.close()


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Process Cigar String..')

    # Required arguments
    parser.add_argument('input_file', help='Input SAM or BAM file.')
    parser.add_argument('IGH_bed_file', help='IGH BED file input.')
    parser.add_argument('IGK_bed_file', help='IGK BED file input.')
    parser.add_argument('IGL_bed_file', help='IGL BED file input.')
    parser.add_argument('species', help='Species name.')

    # Optional arguments
    parser.add_argument('-o', '--output', help='Output directory path.')

    # Parse arguments
    args = parser.parse_args()

    if args.output is None:
        args.output = f"/home1/zhuyixin/sc1/AssmQuality/errorStats/{args.species}/"
        warnings.warn("No output dir specified. Using default output dir path: " + args.output)


    # Validate the input_file extension
    valid_extensions = ['.sam', '.bam']
    file_extension = os.path.splitext(args.input_file)[1]
    if file_extension not in valid_extensions:
        parser.error("Input file must be a SAM or BAM file.")


    if args.input_file.endswith('.bam'):
        process_bam_file(args.input_file, args.IGH_bed_file, args.IGK_bed_file, args.IGL_bed_file, args.output)
        #tree = load_bed_to_interval_tree(args.IGK_bed_file)
        #assert tree['scaffold_15'].overlaps(26968008, 26970008)

if __name__ == "__main__":
    main()
    