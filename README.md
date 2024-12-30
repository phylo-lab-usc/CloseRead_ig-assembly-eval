# IG Assembly Evaluation with CloseRead

**Brief description:** This tool aims to evaluate the assembly of immunoglobulin (IG) sequences. The evaluation process includes various scripts and workflows to analyze, compare, and validate IG sequence assemblies.

![pipeline overview](plots/overview.png)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)

## Installation

Follow these steps to set up the development environment:

```bash
# Clone the repository
git clone https://github.com/phylo-lab-usc/CloseRead_ig-assembly-eval.git

# Navigate to the project directory
cd CloseRead_ig-assembly-eval

# Create and activate the conda environments
conda env create -f ig-assembly-eval.yml
conda activate ig-assembly-eval

# Install closeread
pip install -e .

# Optional: update your system path
```

### Other Requirements

- [IgDetective](https://github.com/Immunotools/IgDetective.git) is required for the following analysis. Please make sure all the relative path is changed to absolute path in IgDetective source code, otherwise error will occur. Also make sure you are running IgDetective with the python in the above conda enviroment. 
- `samtools` is required for the following analysis.
- `python=3.10` is used, should be installed by default in the above conda enviroments. Again make sure whenever running python, you are using the python installed in the above conda enviroment. You can check by runnning `which python`

## Usage
### 1. Running closeread Pipeline

#### Required input files:

- HiFi fastq/BAM files that generated the assembly of the species of interest at `$HOME/$fastqdir/$species_name/`
- Merged diploid assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.merged.fasta`, use the ENTIRE assembly even if you are only looking at specific loci
- Primary/Haplotype1/Maternal assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.pri.fasta`
- Alternate/Haplotype2/Paternal assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.alt.fasta`
- Above assembly files' index file `.fai`
- (Optional, if want to skip IgDetective) Loci Annotation file in `${species_name}.customIG.txt`, file format see example

#### Optional input file:
- `species_metainfo.csv` containing meta information of the species of interest, format see example
- Gene level annotation file in either IgDetective generated format OR in [Gene, Chromosome, Strand, Start, End] csv format

```bash
closeread-pipeline --species species1,species2 \
                   --home /path/to/home \
                   --haploid True \
                   --fastqdir /path/to/fastqdir \
                   --closeread /path/to/closeread \
                   --igdetective_home /path/to/igdetective_home
```

```bash
closeread-plot [OPTIONS]
usage: closeread-plot [-h] (--s species | --sf species_file) --g gene --dirStat errorStatsDir --dirPlot
                    errorPlotsDir [--cov lowCov_threshold] [--p padding] [--re single_read_error]
                    [--rc readview_correct_threshold] [--bc baseview_correct_threshold] [--m meta] [--so stats_only]
                    [--pg gene_level assessment]

CloseRead Evaluation Stats and Visualization.

options:
  -h, --help            show this help message and exit

required:
  --s species           Single species identifier (use this if you are providing one species)
  OR
  --sf species_file     Path to file containing a list of species (use this if you are providing multiple species)

  --g gene              Gene identifier(IGH/IGK/IGL)
  --dirStat errorStatsDir
                        Path to the previous errorStats directory containing mpileup file
  --dirPlot errorPlotsDir
                        Path to the output errorPlots directory
optional:
  --ha haploid          Haploid status (True/False) (alternate IG loci will not be shown if too short or has multiple)
  --m meta              Absolute path to the meta information .csv file, used for generating pdf.
  --so stats_only       output .txt and .csv files only, skip visualization.
  --pg gene_level assessment
                        Absolute path for gene level annotation file, Generate gene level read support information only. 
  --cov lowCov_threshold
                        Threshold for low coverage (default: 2)
  --p padding           Padding around low coverage regions (default: 2000bps)
  --re single_read_error
                        Threshold for a single read to consider as high mismatch (default: 0.01)
  --rc readview_correct_threshold
                        Number of high mismatch reads needed to cover a position for it to be considered as high mismatch
                        rate position from read-view (default: 5)
  --bc baseview_correct_threshold
                        Threshold for the percent of reads with exact match at a position for it to be considered as well-
                        supported, used in heatmap (default: 80 percent)
```

### 2. Output and Intermediate files
#### The intermediate stats files will be in the `errorStats/` directory and should include the following 11 files:

- `IGH.txt`, read-oriented stats for both primary and alternate assembly at IGH locus
- `IGH_alt_pileup.txt`, mpileup file, basepair-oriented stats for alternate assembly at IGH locus
- `IGH_pri_pileup.txt`, mpileup file, basepair-oriented stats for primary assembly at IGH locus
- `IGK.txt`
- `IGK_alt_pileup.txt`
- `IGK_pri_pileup.txt`
- `IGL.txt`
- `IGL_pri_pileup.txt`
- `cigar.end`, empty flag file
- `pileup.end`, empty flag file
- `nonIG.txt`, read-oriented stats for both primary and alternate assembly at non-IG locus

#### The final output files will be in the `errorPlots/` directory and should include the following files:

- `summary.allreads.png`, summary evaluation of the general stats
- `length.png`, locus length plot
- `readcoverage.all.png`, coverage by read mapping quality across a genomic region, for both haplotype if diploid
- `basecoverage.PerCorrect.png`, basepair level coverage (% mismatch per position) and a heatmap of poorly supported positions, for both haplotype if diploid
- `break.txt` containing information for all low coverage positions
- `read.mismatch.txt` containing positional information for all high-mismatch rate region from read-view
- `base.exactmismatch.csv` containing detailed basepair-view mismatch information for every basepair with < bc (80% default) base support
- `base.avgmismatch.csv` containing averaged basepair-view mismatch information, grouping consecutive basepairs (< 1000bps spaced), best for identifying regions with lots of single base mismatches
- `finalMismatch.csv` containing positions marked highly mismatched from both basepair-view and read-view, for both haplotype if diploid
- `result.pdf`, single pdf with all figures and meta information. 
- `{species}.{gene}.genelevel.csv`, gene level quality assessment, if flag `--pg` given

#### For more information on how to interpret the result please refer to this [document](https://docs.google.com/document/d/1QOh3Z6noqZ7x-u70QhQv4VhQOrCJJ1hz_NEmvBDaT_A/pub) 

### 3. Log files
Log files will be avaliable at `$HOME/logs/`

### Folder Structure

A brief example overview of the project's structure and directories:

```plaintext
$HOME/
├── assemblies/          # Input assembly files
│   ├── mEubGla1.merged.fasta       # Example merged assembly
│   ├── mEubGla1.merged.fasta.fai   # Example merged assembly index
│   ├── mEubGla1.pri.fasta          # Example Primary/Haplotype1/Maternal assembly
│   ├── mEubGla1.pri.fasta.fai      # Example Primary/Haplotype1/Maternal assembly index
│   ├── mEubGla1.alt.fasta          # Example Alternate/Haplotype2/Paternal assembly
│   ├── mEubGla1.alt.fasta.fai      # Example Alternate/Haplotype2/Paternal assembly index
├── $fastqdir/           # Input fastq files
│   ├── mEubGla1              # Example dir containing fastq/bam files for mEubGla1
├── errorStats/          # Mismatch and coverage stats files
│   ├── mEubGla1              # Example dir containing output stats files for mEubGla1
├── errorPlots/          # Final visualizations
│   ├── mEubGla1              # Example dir containing output visualization files for mEubGla1
```

## Project Structure

```plaintext
CloseRead_ig-assembly-eval/
├── README.md          # Project documentation
├── ig-assembly-eval.yml       # Conda environment file
├── closeread/              # Source code and scripts
├── plots/             # Directory for plots and figures
├── curated_IGH/       # Directory for LJA curated IGH assembly
├── example/           # Example input format
├── test/           # Test case (TBD)
```

## Citing:

Assessing Assembly Errors in Immunoglobulin Loci: A Comprehensive Evaluation of Long-read Genome Assemblies Across Vertebrates
Yixin Zhu, Corey Watson, Yana Safonova, Matt Pennell, Anton Bankevich
bioRxiv 2024.07.19.604360; doi: https://doi.org/10.1101/2024.07.19.604360

## License

This project is licensed under the GNU General Public License v3 (GPLv3). See the [LICENSE](LICENSE) file for full details.
