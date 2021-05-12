#!/usr/bin/env python3
"""Copy the openff-toolkit examples suite to a local directory"""

from pathlib import Path
from os import environ
from shutil import copytree, ignore_patterns
import argparse


class NoExamplesLibraryError(Exception):
    """Slightly more descriptive error for when the examples library is not found"""

    def __init__(self):
        message = "\n".join(
            (
                "Could not find examples library.",
                "",
                "Is the openff-toolkit-examples Conda package installed?",
            )
        )
        super().__init__(message)


def main():
    """Entrypoint for the script"""
    parser = argparse.ArgumentParser(description=__doc__)

    default_target = "./examples/"
    parser.add_argument(
        "--target",
        default=default_target,
        type=Path,
        help=f"Target directory to store the examples suite. Defaults to {default_target}",
    )

    parser.add_argument(
        "--include-deprecated",
        action="store_true",
        help="Copy over deprecated examples in addition to current examples",
    )

    args = parser.parse_args()

    if args.include_deprecated:
        copytree_ignore = None
    else:
        copytree_ignore = ignore_patterns("deprecated")

    target_path = args.target
    library_path = Path(environ["CONDA_PREFIX"]) / "share/openff-toolkit/examples"

    try:
        copytree(library_path, target_path, ignore=copytree_ignore)
    except FileNotFoundError as err:
        if Path(err.filename) == library_path and err.filename2 is None:
            raise NoExamplesLibraryError() from err
        raise err


if __name__ == "__main__":
    main()
