#!/usr/bin/env python3
"""Copy the openff-toolkit examples suite to a local directory"""

from pathlib import Path
from os import environ
from shutil import copytree


def main():
    """Entrypoint for the script"""
    examples_path = Path(environ["CONDA_PREFIX"]) / "share/openff-toolkit/examples"
    target_path = Path("./examples")

    copytree(examples_path, target_path)


if __name__ == "__main__":
    main()
