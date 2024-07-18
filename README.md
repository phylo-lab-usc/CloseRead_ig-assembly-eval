# IG Assembly Evaluation

**Brief description:** This project aims to evaluate the assembly of immunoglobulin (IG) sequences. The evaluation process includes various scripts and workflows to analyze, compare, and validate IG sequence assemblies.

![pipeline overview](plots/overview.png)

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)

## Installation

Follow these steps to set up the development environment:

```bash
# Clone the repository
git clone https://github.com/yourusername/ig-assembly-eval.git

# Navigate to the project directory
cd ig-assembly-eval

# Create and activate the conda environment
conda env create -f assembly.yml
conda activate ig-assembly-eval
```

### Other Requirements

[IgDetective](https://github.com/Immunotools/IgDetective.git) is required for the following analysis.

## Usage
### Running the Pipeline

![snakefile overview](plots/dag_snake.png)

Use the Snakefile to run all the code located in the `code` folder. Above is an example pipeline overview for 1 species.

#### Required input files:

- HiFi fastq/BAM files of species of interest at `$HOME/$fastqdir/$species_name/`
- Merged diploid assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.merged.fasta`
- Primary/Haplotype1/Maternal assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.pri.fasta`
- Alternate/Haplotype2/Paternal assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.alt.fasta`
- Above assembly files' index file `.fai`

#### Optional input file:
- `species_metainfo.csv` containing meta information of the species of interest

#### Please make sure you modify the header lines of `Snakefile` to reflect your directory organization:

- `SPECIES = ["mEubGla1"]`, list of species name
- `fastqdir = ["hifi_fastq"]`,  sub-directory of your home directory where your fastq files are located
- `HAPLOID = ["False"]`,  if the list of species are halpid or not
- `HOME = "/home1/zhuyixin/zhuyixin_proj/AssmQuality"`,  your home directory

#### The output stats files will be in the `errorStats/` directory and should include the following 11 files:

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

```bash
# Run the main workflow using Snakemake
snakemake -R all --snakefile Snakefile --printshellcmds --reason --verbose --latency-wait 60000 --cores all
```
### Generating Visualizations

After running the Snakefile, execute the evaluation notebook `code/evaluate.ipynb` to generate visualizations, stored in `errorPlots`.

Detailed break analysis is in notebook `code/break.ipynb`.

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
ig-assembly-eval/
├── README.md          # Project documentation
├── Snakefile          # Main workflow file
├── assembly.yml       # Conda environment file
├── code/              # Source code and scripts
├── plots/             # Directory for plots and figures
├── curated_IGH/       # Directory for LJA curated IGH assembly
```
