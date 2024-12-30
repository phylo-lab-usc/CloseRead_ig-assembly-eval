import os
import subprocess
from datetime import datetime

def data_prep(species, home, fastqdir, haploid, closeread):
    """Run data preparation."""
    # Define log file
    log_file = os.path.join(home, "logs", f"data_prep_{species}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define the script path
    script = os.path.join(closeread, "scripts/dataPrepAutomated.sh")

    # Build the command
    cmd = [
        script,
        "-s", species,
        "-w", fastqdir,
        "-h", haploid,
        "-d", home,
    ]

    # Open log file for writing
    with open(log_file, "w") as log:
        # Run the command, redirecting stdout and stderr to the log file
        try:
            subprocess.run(cmd, stdout=log, stderr=log, check=True)
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Data preparation completed successfully. Logs are in {log_file}")
        except subprocess.CalledProcessError as e:
            print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Data preparation failed. Check logs for details: {log_file}")
            raise

if __name__ == "__main__":
    import argparse

    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Run the data preparation step.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--fastqdir", required=True, help="Path to the FASTQ directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True or False).")
    parser.add_argument("--closeread", required=True, help="Path to the CloseRead directory.")

    args = parser.parse_args()

    # Run the data preparation function with parsed arguments
    data_prep(args.species, args.home, args.fastqdir, args.haploid, args.closeread)
