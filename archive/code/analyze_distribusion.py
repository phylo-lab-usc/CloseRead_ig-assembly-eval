import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import seaborn as sns

def analyze_distribution(file_path, threshold, species, chr):
    fractions = np.loadtxt(file_path, usecols=[5])

    # Fit a Gaussian mixture model
    gmm = GaussianMixture(n_components=2)
    gmm.fit(fractions.reshape(-1, 1))

    output_file = f"/home1/zhuyixin/sc2/ImmAssm/cigar_result/weights_{species}_per_threshold.txt"

    # Append the weights to a file
    with open(output_file, "a") as f:
        f.write(f"Threshold {threshold}: {chr}: {gmm.weights_}\n")
    
    fig, (ax_hist, ax_kde) = plt.subplots(1, 2, figsize=(12, 6))

    # Plotting the histogram
    ax_hist.hist(fractions, bins=30, density=True)
    ax_hist.set_title(f'Histogram (Threshold {threshold})')
    ax_hist.set_xlabel('Fraction')
    ax_hist.set_ylabel('Density')

    # Plotting the KDE
    sns.kdeplot(fractions, ax=ax_kde)
    ax_kde.set_title(f'KDE Plot (Threshold {threshold})')
    ax_kde.set_xlabel('Fraction')
    ax_kde.set_ylabel('Density')

    # Saving the histogram
    histogram_output_file = f"/home1/zhuyixin/sc2/ImmAssm/cigar_result/histogram_{species}_{threshold}.png"
    plt.savefig(histogram_output_file)
    plt.close()

    
if __name__ == "__main__":
    analyze_distribution(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])