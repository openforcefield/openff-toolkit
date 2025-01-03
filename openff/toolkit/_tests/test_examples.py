"""
Test that the examples in the repo run without errors.
"""

import pathlib
import re
import subprocess
import tempfile
import textwrap

import pytest

from openff.toolkit._tests.utils import _get_readme_path
from openff.toolkit.utils import RDKIT_AVAILABLE, get_data_file_path, temporary_cd


def run_script_file(file_path):
    """Run through the shell a python script."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        with temporary_cd(tmp_dir):
            cmd = ["python", file_path]
            if "conformer_energies.py" in file_path:
                cmd.append("--filename")
                mol_file = get_data_file_path("molecules/ruxolitinib_conformers.sdf")
                cmd.append(mol_file)
            try:
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError:
                raise Exception(f"Example {file_path} failed")


def run_script_str(script_str):
    """Execute a Python string through the shell in a temporary directory.

    With respect to eval, this has the advantage of catching all import errors.

    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_file_path = (pathlib.Path(tmp_dir) / "temp.py").as_posix()
        # Create temporary python script.
        with open(temp_file_path, "w") as f:
            f.write(script_str)
        # Run the Python script.
        try:
            run_script_file(temp_file_path)
        except Exception as error:
            script_str = textwrap.indent(script_str, "    ")
            raise Exception(f"The following script failed:\n{script_str}") from error


def find_example_scripts() -> list[str]:
    """Find all Python scripts, excluding Jupyter notebooks, in the examples folder.

    Returns
    -------
    example_file_paths : list[str]
        List of full paths to python scripts to execute.
    """
    if "site-packages" in __file__:
        # This test file is being collected from the installed package, which
        # does not provide the examples folder in the same location
        return list()

    examples_dir_path = pathlib.Path(__file__).parents[3] / "examples"

    # Examples that require RDKit
    rdkit_examples = {
        examples_dir_path / "conformer_energies/conformer_energies.py",
    }

    example_file_paths = []
    for example_file_path in examples_dir_path.glob("*/*.py"):
        if not RDKIT_AVAILABLE:
            if example_file_path in rdkit_examples:
                continue
        example_file_paths.append(example_file_path.as_posix())

    return example_file_paths


def find_readme_examples() -> list[str]:
    """Yield the Python scripts in the main README.md file.

    Returns
    -------
    readme_examples : list[str]
        The list of Python scripts included in the README.md files.
    """
    readme_file_path = _get_readme_path()

    if readme_file_path is None:
        return list()

    with open(readme_file_path) as f:
        readme_content = f.read()

    return re.findall("```python(.*?)```", readme_content, flags=re.DOTALL)


@pytest.mark.parametrize("readme_example_str", find_readme_examples())
def test_readme_examples(readme_example_str):
    """Test the example"""
    run_script_str(readme_example_str)


@pytest.mark.parametrize("example_file_path", find_example_scripts())
def test_examples(example_file_path):
    """Test that the example run without errors."""
    run_script_file(example_file_path)
