import subprocess
import os
from datetime import datetime

def convert_primary_bam(species, home):
    """Convert merged BAM to primary BAM."""
    # Define the log file path
    log_file = os.path.join(home, "logs", f"{species}_convert_primary_bam.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define input and output file paths
    bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted.bam")
    output_bam = os.path.join(home, "aligned_bam", species, f"{species}_merged_sorted_primary.bam")
    output_index = f"{output_bam}.csi"

    # Skip step if output already exists
    if os.path.exists(output_bam) and os.path.exists(output_index):
        with open(log_file, "a") as log:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Primary BAM and index already exist: {output_bam}, {output_index}\n")
            log.flush()
        return

    # Run samtools commands and log output
    with open(log_file, "w") as log:
        try:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Converting BAM to primary BAM for species: {species}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Input BAM: {bam}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Output BAM: {output_bam}\n")
            log.flush()
            # Run samtools view
            subprocess.run(
                ["samtools", "view", "-b", "-F", "0x800", "-F", "0x100", "-@", "30", bam, "-o", output_bam],
                stdout=log,
                stderr=log,
                check=True,
            )
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - samtools view completed successfully.\n")
            log.flush()
            # Run samtools index
            subprocess.run(
                ["samtools", "index", "-c", output_bam],
                stdout=log,
                stderr=log,
                check=True,
            )
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - samtools index completed successfully.\n")

        except subprocess.CalledProcessError as e:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during BAM conversion: {e}\n")
            raise RuntimeError(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Failed to convert BAM for species {species}. Check the log: {log_file}") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert merged BAM to primary BAM.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")

    args = parser.parse_args()

    # Call the function
    convert_primary_bam(args.species, args.home)
