import os
import sys
import subprocess
import logging
from .steps.data_prep import data_prep
from .steps.convert_primary_bam import convert_primary_bam
from .steps.loci_location import loci_location
from .steps.final_ig_loci import final_ig_loci
from .steps.cigar_processing import cigar_processing
from .steps.coverage_analysis import coverage_analysis
import argparse

def run_pipeline_cli():
    parser = argparse.ArgumentParser(description="Run the CloseRead pipeline.")
    parser.add_argument("--species", required=True, help="Comma-separated list of species (e.g., species1,species2).")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True or False).")
    parser.add_argument("--fastqdir", required=True, help="Path to the FASTQ directory.")
    parser.add_argument("--closeread", required=True, help="Path to the CloseRead directory.")
    parser.add_argument("--igdetective_home", required=True, help="Path to the igDetective home directory.")
    
    args = parser.parse_args()
    run_pipeline(args)


def configure_logging(home):
    """Configure logging for the pipeline."""
    log_dir = os.path.join(home, "logs")
    os.makedirs(log_dir, exist_ok=True)  # Create the logs directory if it doesn't exist

    logging.basicConfig(
        filename=os.path.join(log_dir, f"pipeline.log"),
        filemode="w",  # Overwrite log file each run
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info(f"Logging initialized. Logs will be written to {os.path.join(log_dir, 'pipeline.log')}")

def run_pipeline(args):
    """Run the pipeline with command-line parameters."""

    # Extract arguments
    species_list = args.species.split(",")
    home = args.home
    haploid = args.haploid
    fastqdir = args.fastqdir
    closeread = args.closeread
    igdetective_home = args.igdetective_home
    configure_logging(home)

    # Validate input directories
    for dir_path, dir_name in [
        (home, "Home directory"),
        (f"{home}/{fastqdir}", "FASTQ directory"),
        (closeread, "CloseRead directory"),
        (igdetective_home, "IGDetective directory"),
    ]:
        if not os.path.exists(dir_path):
            logging.error(f"{dir_name} not found: {dir_path}")
            raise FileNotFoundError(f"{dir_name} not found: {dir_path}")

    logging.info("Starting the pipeline...")

    for species in species_list:
        logging.info(f"Processing species: {species}")
        try:
            # Step 1: Data Preparation
            logging.info(f"Step 1: Data Preparation for {species}")
            data_prep(species, home, fastqdir, haploid, closeread)

            # Step 2: Convert Primary BAM
            logging.info(f"Step 2: Convert Primary BAM for {species}")
            convert_primary_bam(species, home)

            # Step 3: Loci Location
            logging.info(f"Step 3: Loci Location for {species}")
            loci_location(species, home, haploid, igdetective_home)

            # Step 4: Final IG Loci
            logging.info(f"Step 4: Final IG Loci for {species}")
            final_ig_loci(species, home, closeread)

            # Step 5: Coverage Analysis
            logging.info(f"Step 5: Coverage Analysis for {species}")
            coverage_analysis(species, home, closeread)

            # Step 6: CIGAR Processing
            logging.info(f"Step 6: CIGAR Processing for {species}")
            cigar_processing(species, home, closeread)

            logging.info(f"Species {species} processed successfully!")

        except Exception as e:
            logging.error(f"Error processing species {species}: {e}", exc_info=True)
            sys.exit(1)  # Exit with error code

    logging.info("Pipeline completed successfully!")



if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run the pipeline.")
    parser.add_argument("--species", required=True, help="Comma-separated list of species.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True/False).")
    parser.add_argument("--fastqdir", required=True, help="Path to the FASTQ directory.")
    parser.add_argument("--closeread", required=True, help="Path to the CloseRead directory.")
    parser.add_argument("--igdetective_home", required=True, help="Path to the IGDetective directory.")

    args = parser.parse_args()

    # Run the pipeline
    run_pipeline(args)
