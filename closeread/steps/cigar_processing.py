import subprocess
import os
from datetime import datetime

def cigar_processing(species, home, closeread):
    """Process CIGAR data."""
    # Define log file path
    log_file = os.path.join(home, "logs", f"cigar_processing_{species}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define paths for script and input files
    script = os.path.join(closeread, "scripts/cigar.py")
    bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
    annotation = os.path.join(home, "gene_position", f"{species}.final.Ig_loci.txt")
    error_dir = os.path.join(home, "errorStats", species)
    os.makedirs(error_dir, exist_ok=True)

    # Open log file and execute the subprocess
    with open(log_file, "w") as log:
        try:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Starting CIGAR processing for species: {species}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - BAM file: {bam}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Annotation file: {annotation}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - errorStats directory: {error_dir}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running script: python {script} {bam} {annotation} {species} {error_dir}\n")
            log.flush()
            # Run the script
            subprocess.run(
                [
                    "python",
                    script,
                    bam,
                    annotation,
                    species,
                    error_dir,
                ],
                stdout=log,
                stderr=log,
                check=True,
            )
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - CIGAR processing completed successfully for species: {species}\n")

        except subprocess.CalledProcessError as e:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during CIGAR processing: {e}\n")
            raise RuntimeError(f"Failed to process CIGAR data for species {species}. Check the log: {log_file}") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process CIGAR data.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--closeread", required=True, help="Path to the CloseRead directory.")

    args = parser.parse_args()

    # Call the function
    cigar_processing(args.species, args.home, args.closeread)
