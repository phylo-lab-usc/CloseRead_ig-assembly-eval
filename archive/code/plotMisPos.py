import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import matplotlib.ticker as ticker

def plot_filtered_reads(input_file_path, output_file_path):
    # Read the data
    data = pd.read_csv(input_file_path, header=None, sep='\t')
    
    # Filter reads with mismatch rate > 0.02
    # Assuming the mismatch rate is in the 7th column (index 6)
    chr_data = data[data[1] == 'SUPER_21']
    filtered_data = chr_data[chr_data[6] > 0.01]

    # Plot occurrences across genome positions
    plt.figure(figsize=(12, 6))
    num_bins = (150026182 - 148125104) // 100 + 1
    sns.histplot(filtered_data[2], bins=num_bins, color='purple')  # Assuming position is in the third column
    plt.title('Occurrences of Reads with Mismatch Rate > 0.01 Across Genome Positions', fontsize=16, fontweight='bold')
    plt.xlabel('Genome Position', fontsize=14, fontweight='bold')
    plt.ylabel('Number of Reads', fontsize=14, fontweight='bold')
    plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
    #plt.xlim([13912342,15556558]) # rheMac
    #plt.xlim([78085623,78740677]) # bovine
    #plt.xlim([73919,2570154]) # mCanLor1
    #plt.xlim([71533516,72192326]) # mCerEla1
    #plt.xlim([184922455,186976783]) # mCynVol1
    #plt.xlim([56846581,59238841]) # mEleMax1
    #plt.xlim([319143269,322209015]) # mMacEug1
    #plt.xlim([304572038,305242180]) # mMonDom1
    plt.xlim([148125104,150026182]) # mNeoNeb1
    #plt.xlim([105586437,106879844]) # HUMAN
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_file_path, format='png', dpi=300)

def main():
    if len(sys.argv) != 3:
        print("Usage: python cigar_mismatch.py <input_file_path> <output_file_path> ")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    plot_filtered_reads(input_file_path, output_file_path)

if __name__ == "__main__":
    main()
