import subprocess
import os
from datetime import datetime

def coverage_analysis(species, home, closeread, annotation):
    """Run coverage analysis."""
    # Define the log file path
    log_file = os.path.join(home, "logs", f"{species}_coverage_analysis.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define input paths and directories
    script = os.path.join(closeread, "scripts/coverage.sh")
    assembly = os.path.join(home, "assemblies", f"{species}.merged.fasta")
    bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
    error_dir = os.path.join(home, "errorStats", species)
    os.makedirs(error_dir, exist_ok=True)

    with open(log_file, "w") as log:
        try:
            # Log the start of the process
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Starting coverage analysis for species: {species}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - BAM file: {bam}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Annotation file: {annotation}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - errorStats directory: {error_dir}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running script: {script} -s {species} -a {assembly} -b {bam} -f {annotation} -d {home}\n")
            log.flush()
            # Run the script
            subprocess.run(
                f"{script} -s {species} -a {assembly} -b {bam} -f {annotation} -d {home}",
                stdout=log,
                stderr=log,
                shell=True,
                check=True,
            )
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Coverage analysis completed successfully for species: {species}\n")
            log.flush()
        except subprocess.CalledProcessError as e:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during coverage analysis: {e}\n")
            raise RuntimeError(f"Failed to complete coverage analysis for species {species}. Check the log: {log_file}") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run coverage analysis.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--closeread", required=True, help="Path to the CloseRead directory.")
    parser.add_argument("--annotation", required=True, help="Path to the annotation file.")

    args = parser.parse_args()

    # Call the function
    coverage_analysis(args.species, args.home, args.closeread, args.annotation)
