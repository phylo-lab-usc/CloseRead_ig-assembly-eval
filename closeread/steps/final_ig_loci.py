import subprocess
import os
from datetime import datetime

def final_ig_loci(species, home, closeread):
    """Process loci into final IG loci."""
    # Define the log file path
    log_file = os.path.join(home, "logs", f"final_ig_loci_{species}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define script and output paths
    script = os.path.join(closeread, "scripts/finalGene.py")
    output = os.path.join(home, "gene_position", f"{species}.final.Ig_loci.txt")

    with open(log_file, "w") as log:
        try:
            # Log start of the process
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Starting final IG loci processing for species: {species}\n")
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Output file: {output}\n")

            # Check if output exists
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running script: python {script} {species} {home}\n")
            log.flush()
            subprocess.run(
                f"python {script} {species} {home}",
                stdout=log,
                stderr=log,
                check=True,
            )
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Final IG loci generated for {species}.\n")


        except subprocess.CalledProcessError as e:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during IG loci processing: {e}\n")
            raise RuntimeError(f"Failed to process final IG loci for species {species}. Check the log: {log_file}") from e


if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process loci into final IG loci.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--closeread", required=True, help="Path to the CloseRead directory.")

    args = parser.parse_args()

    # Call the function
    final_ig_loci(args.species, args.home, args.closeread)