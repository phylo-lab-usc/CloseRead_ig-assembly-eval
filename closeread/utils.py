import yaml
import os
import logging

def setup_logging(log_file):
    """
    Set up logging for a step.
    Args:
        log_file (str): Path to the log file.
    """
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        filemode="w",
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logging.info(f"Logging initialized. Writing to {log_file}")


def extract_requirements(yaml_file, output_file):
    with open(yaml_file, "r") as f:
        env = yaml.safe_load(f)

    dependencies = env.get("dependencies", [])
    pip_dependencies = []
    conda_dependencies = []

    for dep in dependencies:
        if isinstance(dep, str):  # Conda dependencies
            conda_dependencies.append(dep)
        elif isinstance(dep, dict) and "pip" in dep:  # Pip dependencies
            pip_dependencies.extend(dep["pip"])

    # Write to requirements.txt
    with open(output_file, "w") as f:
        for dep in conda_dependencies:
            if "=" in dep:  # Keep Python libraries only
                if not any(dep.startswith(non_py) for non_py in ["_", "perl", "r-", "lib", "xorg", "sqlite"]):
                    f.write(dep.replace("=", "==") + "\n")
        for dep in pip_dependencies:
            f.write(dep + "\n")

    print(f"requirements.txt created at {output_file}")

# Run the script
extract_requirements("/home1/zhuyixin/zhuyixin_proj/AssmQuality/code/CloseRead_ig-assembly-eval/ig-assembly-eval.yml", "/home1/zhuyixin/zhuyixin_proj/AssmQuality/code/CloseRead_ig-assembly-eval/requirements.txt")
