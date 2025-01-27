# IG Assembly Evaluation with CloseRead

**Brief description:** This tool aims to evaluate the assembly of immunoglobulin (IG) sequences. The evaluation process includes various scripts and workflows to analyze, compare, and validate IG sequence assemblies.

![pipeline overview](plots/overview.png)

## Table of Contents

- [Installation](#installation)
- [Test Case](#test-case)
- [Input](#input)
- [Running](#running)
- [Output and Intermediate files](#output)
- [Folder Structure](#folder-structure)
- [Project Structure](#project-structure)
- [Citing](#citing)

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
- `samtools` is required for the following analysis
- `python=3.10` is used, should be installed by default in the above conda enviroments. Again make sure whenever running python, you are using the python installed in the above conda enviroment. You can check by runnning `which python`

## Test Case
To check if CloseRead was correctly intalled, we suggest running the following test case prior to any of your analysis. Run:
```
closeread-pipeline --species test --home {PATH to Closeread}/test/ --haploid True --fastqdir hifi_fastq --closeread /{PATH to Closeread}/closeread --igdetective_home {PATH to IgDetective}
closeread-plot --s test --g IGH --home {PATH to Closeread}/test/ --dirStat {PATH to Closeread}/test/errorStats --dirPlot {PATH to Closeread}/test/errorPlots/ 
```
If using 32 cores, the entire process should be finished within 5 minutes. 
The output can be found in `{PATH to Closeread}/test/`, we also provide the expected output here [link](https://figshare.com/s/9b8db110af1511871669). Please compare the figures and any outputs to make sure they are the same.
For more information on how to interpret the result please refer to this [document](https://docs.google.com/document/d/1QOh3Z6noqZ7x-u70QhQv4VhQOrCJJ1hz_NEmvBDaT_A/pub)


## Input
#### Required input files:

- HiFi fastq/BAM files that generated the assembly of the species of interest at `$HOME/$fastqdir/$species_name/`
- Merged diploid assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.merged.fasta`, use the ENTIRE assembly even if you are only looking at specific loci
- Primary/Haplotype1/Maternal assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.pri.fasta`
- Alternate/Haplotype2/Paternal assembly fasta file of species of interest at `$HOME/assemblies/${species_name}.alt.fasta`
- Above assembly files' index file `.fai`
- (Optional, if loci position already known and want to skip IgDetective) Loci Annotation file in `${species_name}.customIG.txt`, **file format see example/Emax.customIG.txt**

#### Optional input file:
- `species_metainfo.csv` containing meta information of the species of interest, format see example
- `Gene level annotation file` in either IgDetective generated format **OR** in [Gene, Chromosome, Strand, Start, End] csv format

## Running
#### Run read-to-assembly pipeline
```bash
closeread-pipeline [OPTIONS]
usage: closeread-pipeline [-h] --species SPECIES --home HOME --haploid HAPLOID --fastqdir FASTQDIR --closeread CLOSEREAD
                          (--igdetective_home IGDETECTIVE_HOME | --customIG CUSTOMIG) [--t T]

options:
  -h, --help            show this help message and exit
  --species SPECIES     Comma-separated list of species (e.g., species1,species2).
  --home HOME           Path to the home directory.
  --haploid HAPLOID     Haploid status (True or False).
  --fastqdir FASTQDIR   Path to the FASTQ directory.
  --closeread CLOSEREAD
                        Path to the CloseRead directory.
  --igdetective_home IGDETECTIVE_HOME
                        Path to the IGDetective directory.
  --t T                 # of threads to use (default: 32).
  --customIG CUSTOMIG   Path to DIRECTORY containing ${species_name}.customIG.txt (not path directly to the file), you could put it in the gene_position folder
```
#### Run evalutaion and plotting pipeline
To note, you will need to run each loci separatly for plotting.
Gene level assessment will also need to be ran separatly with `--pg` parameter provided.
```bash
closeread-plot [OPTIONS]
usage: closeread-plot [-h] (--s species | --sf species_file) --g gene --home home --dirStat errorStatsDir --dirPlot errorPlotsDir [--cov lowCov_threshold]
                      [--p padding] [--re single_read_error] [--rc readview_correct_threshold] [--bc baseview_correct_threshold] [--m meta]
                      [--so stats_only] [--pg gene_level assessment]

Run CloseRead Final Evalutation and Plotting.

options:
  -h, --help            show this help message and exit
  --s species           Single species identifier (use this if you are providing one species)
  --sf species_file     Path to file containing a list of species (use this if you are providing multiple species)
  --g gene              Gene identifier(IGH/IGK/IGL)
  --home home           Path to the home directory
  --dirStat errorStatsDir
                        Path to the previous errorStats directory containing mpileup file
  --dirPlot errorPlotsDir
                        Path to the output errorPlots directory
  --cov lowCov_threshold
                        Threshold for low coverage (default: 2)
  --p padding           Padding around low coverage regions (default: 2000bps)
  --re single_read_error
                        Threshold for a single read to consider as high mismatch (default: 0.01)
  --rc readview_correct_threshold
                        Number of high mismatch reads needed to cover a position for it to be considered as high mismatch rate position from read-view
                        (default: 5)
  --bc baseview_correct_threshold
                        Threshold for the percent of reads with exact match at a position for it to be considered as well-supported, used in heatmap
                        (default: 80 percent)
  --m meta              Absolute path to the meta information .csv file, used for generating pdf.
  --so stats_only       output .txt and .csv files only, skip visualization.
  --pg gene_level assessment
                        Absolute path for gene level annotation file, Generate gene level read support information only.
```

## Output and Intermediate files
#### `closeread-pipeline` will generate intermediate stats files in the `{HOME}/errorStats/` directory and should include the following 11 files:

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

#### Other intermediate files generated by `closeread-pipeline`:

- `{HOME}/aligned_bam/{species}/` will contain the aligned bam file and the filtered (only primary alignment) bam file, could be used for IGV visualizations
- `{HOME}/aligned_sam/{species}/` will contain the aligned sam file
- `{HOME}/igGene/{species}.{pri/alt}.igdetective/` will contain the IgDetective output, gene position information can be found here
- `{HOME}/gene_position/{species}.final.Ig_loci.txt` will contain the formated IG loci positions found by IgDetective. If custom IG loci used, for simplicity, custom IG loci files can also be put here and you just need to provide `--customIG {HOME}/gene_position/`

#### `closeread-plot` will use above files as input and generate the final output files in the `{HOME}/errorPlots/` directory and should include the following files:

- `summary.allreads.png`, summary evaluation of the general stats
- `length.png`, locus length plot
- `readcoverage.all.png`, coverage by read mapping quality across a genomic region, for both haplotype if diploid
- `basecoverage.PerCorrect.png`, basepair level coverage (% mismatch per position) and a heatmap of poorly supported positions, for both haplotype if diploid
- `break.txt` containing information for all low coverage positions
- `read.mismatch.txt` containing positional information for all high-mismatch rate region from read-view
- `base.exactmismatch.csv` containing detailed basepair-view mismatch information for every basepair with < bc (80% default) base support
- `base.avgmismatch.csv` containing averaged basepair-view mismatch information, grouping consecutive basepairs (< 1000bps spaced), best for identifying regions with lots of single base mismatches
- `finalMismatch.csv` containing positions marked highly mismatched from both basepair-view and read-view, for both haplotype if diploid
- `result.pdf`, single pdf with all figures and meta information, if meta info given
- `{species}.{gene}.genelevel.csv`, gene level quality assessment, if flag `--pg` given

#### For more information on how to interpret the result please refer to this [document](https://docs.google.com/document/d/1QOh3Z6noqZ7x-u70QhQv4VhQOrCJJ1hz_NEmvBDaT_A/pub) 

#### Log files
Log files will be avaliable at `$HOME/logs/`


## Folder Structure

A brief example overview of the working directorie's structure:

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
├── aligned_sam/          # Aligned sam file
│   ├── mEubGla1       
├── aligned_bam/          # Aligned bam files
│   ├── mEubGla1       
│   │   ├── mEubGla1_merged_sorted.bam          # Aligned bam file
│   │   ├── mEubGla1_merged_sorted.bam.csi          # Aligned bam file index
│   │   ├── mEubGla1_merged_sorted_primary.bam          # Filtered (only primary alignment) bam file
│   │   ├── mEubGla1_merged_sorted_primary.bam.csi          # Filtered (only primary alignment) bam file index
├── igGene/          # IgDetective Result
│   ├── mEubGla1.pri.igdetective/    
│   ├── mEubGla1.alt.igdetective/          # If diploid 
├── gene_position/          # IG loci position
│   ├── mEubGla1.final.Ig_loci.txt or mEubGla1.customIG.txt   
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
├── test/           # Test case
```

## Citing:

Assessing Assembly Errors in Immunoglobulin Loci: A Comprehensive Evaluation of Long-read Genome Assemblies Across Vertebrates
Yixin Zhu, Corey Watson, Yana Safonova, Matt Pennell, Anton Bankevich
bioRxiv 2024.07.19.604360; doi: https://doi.org/10.1101/2024.07.19.604360

## License

This project is licensed under the GNU General Public License v3 (GPLv3). See the [LICENSE](LICENSE) file for full details.

## Contributing

We loved to hear feedback on this project! Before starting, please check out our [CONTRIBUTING.md](./CONTRIBUTING.md) file for guidelines on how to help improve CloseRead effectively.

**Note:** We prefer issues to be submitted first before creating pull requests. This ensures alignment with the project goals and avoids duplicate efforts.
