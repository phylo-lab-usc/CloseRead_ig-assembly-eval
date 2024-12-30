import subprocess
import os
from datetime import datetime

def loci_location(species, home, haploid, igdetective_home):
    """Run loci location detection for primary and alternate genomes."""
    # Define log file
    log_file = os.path.join(home, "logs", f"loci_location_{species}.log")
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    # Define input and output paths
    pri_genome = os.path.join(home, "assemblies", f"{species}.pri.fasta")
    alt_genome = os.path.join(home, "assemblies", f"{species}.alt.fasta")
    pri_outdir = os.path.join(home, "igGene", f"{species}.pri.igdetective")
    alt_outdir = os.path.join(home, "igGene", f"{species}.alt.igdetective")

    with open(log_file, "w") as log:
        try:
            # Run for primary genome
            if not os.path.exists(pri_outdir):
                os.makedirs(pri_outdir, exist_ok=True)
                log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running loci detection for primary genome: {pri_genome}\n")
                log.flush()
                subprocess.run(
                    [
                        "python",
                        os.path.join(igdetective_home, "run_iterative_igdetective.py"),
                        pri_genome,
                        pri_outdir,
                    ],
                    stdout=log,
                    stderr=log,
                    check=True,
                )
                log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Primary loci detection completed successfully. Output: {pri_outdir}\n")
                log.flush()
            else:
                log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Primary loci output already exists: {pri_outdir}\n")
                log.flush()

            # Run for alternate genome if haploid is False
            if haploid == "False" and not os.path.exists(alt_outdir):
                os.makedirs(alt_outdir, exist_ok=True)
                log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Running loci detection for alternate genome: {alt_genome}\n")
                log.flush()
                subprocess.run(
                    [
                        "python",
                        os.path.join(igdetective_home, "run_iterative_igdetective.py"),
                        alt_genome,
                        alt_outdir,
                    ],
                    stdout=log,
                    stderr=log,
                    check=True,
                )
                log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Alternate loci detection completed successfully. Output: {alt_outdir}\n")
                log.flush()
            elif haploid == "False":
                log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Alternate loci output already exists: {alt_outdir}\n")
                log.flush()

        except subprocess.CalledProcessError as e:
            log.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Error during loci detection: {e}\n")
            raise RuntimeError(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - Loci detection failed for species {species}. Check the log: {log_file}") from e



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run loci location detection.")
    parser.add_argument("--species", required=True, help="Species name.")
    parser.add_argument("--home", required=True, help="Path to the home directory.")
    parser.add_argument("--haploid", required=True, help="Haploid status (True or False).")
    parser.add_argument("--igdetective_home", required=True, help="Path to the igDetective home directory.")

    args = parser.parse_args()

    loci_location(args.species, args.home, args.haploid, args.igdetective_home)
