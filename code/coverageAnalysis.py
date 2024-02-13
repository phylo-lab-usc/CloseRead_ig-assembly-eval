import re
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
import pandas as pd


def parse_read_bases(read_bases):
    """
    Parse read bases to count correct matches (., and ,) while correctly handling indels,
    read starts (^), read ends ($), and skipping over sequences indicating indels (+nXXX or -nXXX).
    """
    # Remove read start markers (^) along with the following character indicating mapping quality
    cleaned_bases = re.sub(r'\^.', '', read_bases)

    # Remove read end markers ($)
    cleaned_bases = cleaned_bases.replace('$', '')

    # Correctly handle indels: remove sequences following + or - indicating the length and actual indel
    cleaned_bases = re.sub(r'[\+\-](\d+)([ACGTNacgtn]+)', '', cleaned_bases)

    # Correctly handle * symbols representing deletions of the reference base (not counted as correct)
    cleaned_bases = cleaned_bases.replace('*', '')

    # After removing all special cases, count '.' and ',' as correct matches
    correct = cleaned_bases.count('.') + cleaned_bases.count(',')

    return correct

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='Analyze per basepair coverage..')

    # Required arguments
    parser.add_argument('input_file', help='Input pileup file.')
    parser.add_argument('output_file', help='Output txt file.')
    parser.add_argument('species', help='Species name.')
    parser.add_argument('gene', help='Gene name.')

    # Parse arguments
    args = parser.parse_args()

    results = {}  # To store results, keyed by (chromosome, position)

    with open(args.input_file, 'r') as f:
        for line in f:
            chrom, pos, ref_base, depth, read_bases, _ = line.split()[:6]
            correct = parse_read_bases(read_bases)
            
            # Calculate percentage of correct calls if depth is not zero to avoid division by zero
            if int(depth) > 0:
                percent_correct = (correct / int(depth)) * 100
            else:
                percent_correct = 0
            
            # Store or process results
            results[(chrom, pos)] = (correct, percent_correct, depth)

    # Output the results to a new file
    with open(args.output_file, 'w') as f:
        for key, (correct, percent_correct, depth) in results.items():
            chrom, pos = key
            f.write(f"{chrom}\t{pos}\t{correct}\t{percent_correct:.2f}%\t{depth}\n")

    print(f"Output written to {args.output_file}")


    window_size = 100  # Define the window size (100 bps)

    # Assuming `results` is your dictionary
    # Convert `results` dictionary to a DataFrame
    data_list = [(*key, *value) for key, value in results.items()]
    data = pd.DataFrame(data_list, columns=['Chrom', 'Pos', 'Correct', 'PercentCorrect', "Depth"])
    data['Pos'] = pd.to_numeric(data['Pos'])
    data.sort_values(by=['Chrom', 'Pos'], inplace=True)

    # Filtering for the first chromosome as an example
    first_chrom = data['Chrom'].unique()[0]
    data_filtered = data[data['Chrom'] == first_chrom]
    data_filtered['RollingAvgCorrect'] = data_filtered['Correct'].rolling(window=window_size, min_periods=1).mean()
    data_filtered['RollingAvgPercent'] = data_filtered['PercentCorrect'].rolling(window=window_size, min_periods=1).mean()
    data_filtered['RollingAvgDepth'] = data_filtered['Depth'].rolling(window=window_size, min_periods=1).mean()

    # Plotting Rolling Average of Correct Matches
    plt.figure(figsize=(10, 6))
    plt.plot(data_filtered['Pos'], data_filtered['RollingAvgCorrect'], label='Rolling Avg Correct Matches', color='tab:green')
    plt.xlabel('Position')
    plt.ylabel('Rolling Avg Correct Matches')
    plt.title(f'Rolling Average of Correct Matches over 100bps for {first_chrom}')
    plt.legend()
    plt.savefig(f'/home1/zhuyixin/sc1/AssmQuality/errorstats/{args.species}/{args.gene}_rolling_avg_correct_matches.png', dpi=300)
    plt.show()

    # Plotting Rolling Average of Percentage Correct
    plt.figure(figsize=(10, 6))
    plt.plot(data_filtered['Pos'], data_filtered['RollingAvgPercent'], label='Rolling Avg Percentage Correct', color='tab:purple')
    plt.xlabel('Position')
    plt.ylabel('Rolling Avg Percentage Correct')
    plt.title(f'Rolling Average of Percentage Correct over 100bps for {first_chrom}')
    plt.legend()
    plt.savefig(f'/home1/zhuyixin/sc1/AssmQuality/errorstats/{args.species}/{args.gene}_rolling_avg_percentage_correct.png', dpi=300)
    plt.show()

    # Plotting Rolling Average of RollingAvgDepth
    plt.figure(figsize=(10, 6))
    plt.plot(data_filtered['Pos'], data_filtered['RollingAvgDepth'], label='Rolling Avg Depth', color='tab:green')
    plt.xlabel('Position')
    plt.ylabel('Rolling Avg Depth')
    plt.title(f'Rolling Average of Depth over 100bps for {first_chrom}')
    plt.legend()
    plt.savefig(f'/home1/zhuyixin/sc1/AssmQuality/errorstats/{args.species}/{args.gene}_rolling_avg_Depth.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()