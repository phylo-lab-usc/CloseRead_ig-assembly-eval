import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import warnings
import pandas as pd
from fpdf import FPDF
import cairosvg
from IPython.display import SVG, display
import math
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from fpdf import FPDF
from PIL import Image
from intervaltree import Interval, IntervalTree
# Import necessary functions from CloseRead_plot.py
from CloseRead_plot import plot_locus_length, plot_summary, plot_coverage, plot_mismatch_coverage

# Import necessary functions from CloseRead_pdf.py
from CloseRead_pdf import make_pdf, make_pdfdi

# Import necessary functions from CloseRead_Stats.py
from CloseRead_stats import read_species_from_file, create_directories
from CloseRead_stats import process_pileup, process_low_coverage_regions, process_read_file, coverage
from CloseRead_stats import calculate_bin_counts, write_pileup, find_overlapping_mismatch_regions


warnings.filterwarnings('ignore')

# Main block to execute
if __name__ == "__main__":
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="CloseRead Evaluation Stats and Visualization.")

    # Add mutually exclusive group for either species or species file input
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--species', type=str, metavar="s", help="Single species identifier (use this if you are providing one species)")
    group.add_argument('--species_file', type=str, metavar="sf", help="Path to file containing a list of species (use this if you are providing multiple species)")

    # Add arguments for gene and haploid
    parser.add_argument('--gene', type=str, required=True, metavar="g", help="Gene identifier(IGH/IGK/IGL)")
    parser.add_argument('--haploid', type=bool, default=False, metavar="h", help="Haploid status (True/False) (alternate IG loci will not be shown if too short or has multiple)")

    # Add argument for Output
    parser.add_argument('--errorStatsDir', type=str, required=True, metavar="dirStat", help="Path to the previous errorStats directory containing mpileup file")
    parser.add_argument('--errorPlotsDir', type=str, required=True, metavar="dirPlot", help="Path to the output errorPlots directory")

    # Add optional argument
    parser.add_argument('--lowCov_threshold', type=int, default=2, metavar="cov", help="Threshold for low coverage (default: 2)")
    parser.add_argument('--padding', type=int, default=2000, metavar="p", help="Padding around low coverage regions (default: 2000bps)")
    parser.add_argument('--single_read_error', type=float, default=0.01, metavar="re", help="Threshold for a single read to consider as high mismatch (default: 0.01)")
    parser.add_argument('--readview_correct_threshold', type=int, default=5, metavar="rc", help="Number of high mismatch reads needed to cover a position for it to be considered as high mismatch rate position from read-view (default: 5)")
    parser.add_argument('--baseview_correct_threshold', type=int, default=5, metavar="bc", help="Threshold for the percent of reads with exact match at a position for it to be considered as well-supported, used in heatmap (default: 80 percent)")
    parser.add_argument('--meta', type=str, metavar="m", help="Absolute path to the meta information .csv file, used for generating pdf.")


    # Parse the command line arguments
    args = parser.parse_args()

    # Determine if single species or species file is provided
    if args.species_file:
        species_list = read_species_from_file(args.species_file)
    else:
        species_list = [args.species]  # Convert single species into a list
   
    gene = args.gene
    haploid = False
    single_read_error = args.single_read_error
    readview_correct_threshold = args.readview_correct_threshold
    chr1_color = "#6AABD7"
    chr2_color = "#F0DDB8"

    # Iterate through each species
    for species in species_list:
        haploid = False
        dirOut = f"{args.errorPlotsDir}/{species}"
        dirStat = f"{args.errorStatsDir}/{species}"
        create_directories(species, dirStat, dirOut)

        #process basepair-view mpileup file
        pri_pileup_file = f'{dirStat}/{gene}_pri_pileup.txt'
        alt_pileup_file = f'{dirStat}/{gene}_alt_pileup.txt'
        pri_pileup = process_pileup(pri_pileup_file)
        # Check if alternate pileup file exists
        if os.path.exists(alt_pileup_file):
            alt_pileup = process_pileup(alt_pileup_file)
            haploid = False
        else:
            print(f"Haploit = True : Alternate mpileup file not found: {alt_pileup_file}")
            alt_pileup = None
            haploid = True
            
        # Concatenate primary and alternate pileup data if diploid
        if alt_pileup is not None:
            merged_pileup = pd.concat([pri_pileup, alt_pileup])
        else:
            merged_pileup = pri_pileup
        merged_pileup.sort_values(by=['Chrom', 'Pos'], inplace=True)

        # Process read-oriented input
        read = process_read_file(f"{dirStat}/{gene}.txt", dirStat)
        chr1 = pri_pileup['Chrom'].value_counts().index[0]
        read_pri = read[(read['chromosome'] == chr1)]
        pri_pileup = pri_pileup[pri_pileup['Chrom']==chr1]
        if not haploid:
            # skipping alternate IG if more than 2 present
            if len(alt_pileup['Chrom'].value_counts()) > 2:
                haploid = True
                chr2 = None
                read_alt = None
            else:
                haploid = False
                chr2 = alt_pileup['Chrom'].value_counts().index[0]
                read_alt = read[(read['chromosome'] == chr2)]
                alt_pileup = alt_pileup[alt_pileup['Chrom']==chr2]

        # Find and output break locations 
        break_regions_list, break_regions = process_low_coverage_regions(merged_pileup, lowCov_threshold=args.lowCov_threshold, padding=args.padding)
        with open(f"{dirOut}/{gene}.break.txt", 'w') as file:
            for row in break_regions_list:
                s = str(row)
                file.write(s+'\n')
        start_break_pri= break_regions[break_regions['Chrom'] == chr1]['Start'].tolist()
        end_break_pri = break_regions[break_regions['Chrom'] == chr1]['End'].tolist()
        if not haploid:
            start_break_alt= break_regions[break_regions['Chrom'] == chr2]['Start'].tolist()
            end_break_alt = break_regions[break_regions['Chrom'] == chr2]['End'].tolist()
            
        #Plot loci length
        plot_locus_length(pri_pileup, alt_pileup, gene, "#6AABD7", "#F0DDB8", dirOut, haploid, chr1, chr2)
        #Plot summary read info
        plot_summary(read_pri, read_alt, "#6AABD7", "#F0DDB8", haploid, dirOut, gene, chr1, chr2)

        #Compute read-view stats
        start_indices_pri, end_indices_pri, high_mismatch_bool_pri, positions_pri, coverage_counts_pri, zero_counts_pri, min_position_pri, max_position_pri, mid_counts_pri, high_mismatch_coverage_pri = coverage(read_pri, single_read_error, readview_correct_threshold)
        min_position_pri = min(pri_pileup['Pos'])
        max_position_pri = max(pri_pileup['Pos'])
        if not haploid:
            start_indices_alt, end_indices_alt, high_mismatch_bool_alt, positions_alt, coverage_counts_alt, zero_counts_alt, min_position_alt, max_position_alt, mid_counts_alt, high_mismatch_coverage_alt = coverage(read_alt, single_read_error, readview_correct_threshold)
            min_position_alt = min(alt_pileup['Pos'])
            max_position_alt = max(alt_pileup['Pos'])

        #Write read-view mismatch to file
        with open(f"{dirOut}/{gene}.read.mismatch.txt", 'w') as file:
            for start, end in zip(start_indices_pri, end_indices_pri):
                file.write(f"{chr1}:{start}-{end}\n")
            if not haploid:
                for start, end in zip(start_indices_alt, end_indices_alt):
                    file.write(f"{chr2}:{start}-{end}\n")

        #Plot read coverage across loci
        plot_coverage(positions_pri, coverage_counts_pri, mid_counts_pri, zero_counts_pri, start_indices_pri,
                      end_indices_pri, high_mismatch_bool_pri.size, start_break_pri, end_break_pri, min_position_pri, max_position_pri, chr1, gene, dirOut)
        if not haploid:
            plot_coverage(positions_alt, coverage_counts_alt, mid_counts_alt, zero_counts_alt, start_indices_alt, 
                          end_indices_alt, high_mismatch_bool_alt.size, start_break_alt, end_break_alt, min_position_alt, max_position_alt, chr2, gene, dirOut)

        #Compute basepair mismatch
        bin_count = calculate_bin_counts(pri_pileup, baseview_correct_threshold=80, bin_size=1000)
        if not haploid:
            alt_bin_count = calculate_bin_counts(alt_pileup, baseview_correct_threshold=80, bin_size=1000)

        #Write base-view mismatch exact/average rough position to file - removes the content if the files already exist
        with open(f"{dirOut}/{gene}.base.exactmismatch.csv", 'w') as f:
            f.write("Chrom,Pos,Correct,PercentCorrect,Depth,Indel\n")
        with open(f"{dirOut}/{gene}.base.avgmismatch.csv", 'w') as f:
            f.write("Chrom,Start,End,AvgPercentMismatch,AvgDepth,AvgIndel\n")
        grouped_baseMis_pri = write_pileup(pri_pileup, gene, dirOut)
        grouped_baseMis_alt = write_pileup(alt_pileup, gene, dirOut)

        #Compute overlapping mismatch regions between read-view and base-view
        overlaps_pri = find_overlapping_mismatch_regions(pri_pileup, start_indices_pri, end_indices_pri)
        overlaps_alt = find_overlapping_mismatch_regions(alt_pileup, start_indices_alt, end_indices_alt)
        overlaps_pri.to_csv(f"{dirOut}/{gene}.{chr1}.finalMismatch.csv", index=False)
        overlaps_alt.to_csv(f"{dirOut}/{gene}.{chr2}.finalMismatch.csv", index=False)

        #Plot basepair mismatch across loci
        #Define the colors for the heatmap
        colors = [(1, 1, 1), (0.6, 0.6, 0.6), (0.2, 0.2, 0.2), (0, 0, 0)]  
        n_bins = 100  
        cmap_name = 'custom'
        cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
        plot_mismatch_coverage(pri_pileup, bin_count, positions_pri, start_indices_pri, end_indices_pri, 
                               high_mismatch_bool_pri.size, start_break_pri, end_break_pri, 
                               chr1_color, chr1, gene, dirOut, cm)
        if not haploid:
            plot_mismatch_coverage(alt_pileup, bin_count, positions_alt, start_indices_alt, end_indices_alt, 
                               high_mismatch_bool_alt.size, start_break_alt, end_break_alt, 
                               chr2_color, chr2, gene, dirOut, cm)

        #Generate PDF - without dotplot
        if args.meta:
            meta = pd.read_csv(args.meta, sep=",")
            LatinName = meta[meta['IndividualID'] == species]['LatinName'].item()
            CommonName = meta[meta['IndividualID'] == species]['CommonName'].item()
            Source = meta[meta['IndividualID'] == species]['Source'].item()
            SourceLink = meta[meta['IndividualID'] == species]['SourceLink'].item()
            Haplotype = meta[meta['IndividualID'] == species]['Haplotype Resolved'].item()
            if Haplotype == "No":
                hapkind = "Not Haplotype Resolved"
            else:
                hapkind = "Haplotype Resolved"
                
            output_filename = f"{dirOut}/{species}_{gene}_result.pdf"
            if not haploid:
                image_files = [
                    f'{dirOut}/{gene}.summary.allreads.png', 
                    f'{dirOut}/{gene}.length.png', 
                    f'{dirOut}/{gene}.{chr1}.readcoverage.all.png',
                    f'{dirOut}/{gene}.{chr1}.basecoverage.PerCorrect.png',
                    f'{dirOut}/{gene}.{chr2}.readcoverage.all.png',
                    f'{dirOut}/{gene}.{chr2}.basecoverage.PerCorrect.png',
                ]
                make_pdfdi(image_files, output_filename, species, CommonName, LatinName, hapkind, Source, overlapx=0.7, overlapy=0.01, scale_top=0.6, scale_bottom=0.3)
            else:
                image_files = [
                    f'{dirOut}/{gene}.summary.allreads.png', 
                    f'{dirOut}/{gene}.length.png', 
                    f'{dirOut}/{gene}.{chr1}.readcoverage.all.png',
                    f'{dirOut}/{gene}.{chr1}.basecoverage.PerCorrect.png',
                ]
                make_pdf(image_files, output_filename, species, CommonName, LatinName, hapkind, Source, overlapx=0.7, overlapy=0.01, scale_top=0.6, scale_bottom=0.3)
                

